
##======================================================================================
## CCISS Publication Scripts
## Step 4h - Figure - Manuscript and GIF plots of persistence vs expansion
##======================================================================================

# Colin Mahony
# c_mahony@alumni.ubc.ca
# 778-288-4008
# July 21, 2019


require (RGtk2)
require(plyr)
require (rChoiceDialogs)
require (data.table)
require(doBy)
require (utils)
require(labdsv)
require(tools )
require(svDialogs)
require(tcltk)
require(randomForest)
require(foreach)
require(dplyr)
require(reshape2)
require(reshape)
library(doParallel)
require(data.table)
library(MASS)   
library(scales)
library(stats)
library(rgl)
library(RColorBrewer)
library(FNN)
library(igraph)
library(raster)
library(maps)
library(mapdata)
library(maptools)
library(sp)
library(colorRamps)
library(rgeos)
library(rgdal)
library(foreign)


#===============================================================================
# Set analysis Parameters
#===============================================================================


setwd("C:\\Colin\\Projects\\2019_CCISS")

grid <- "BC2kmGrid"
hist.years <- c(1995, 2004, 2005, 2009, 2014, 2017)
hist.year.name <- c("1991-2000", "1991-2017", "2001-2010", "2001-2017","2011-2017", "2017")

GCMs <-  c("ACCESS1-0","CanESM2","CCSM4","CESM1-CAM5","CNRM-CM5","CSIRO-Mk3-6-0", "GFDL-CM3","GISS-E2R", "HadGEM2-ES", "INM-CM4", "IPSL-CM5A-MR", "MIROC-ESM", "MIROC5", "MPI-ESM-LR","MRI-CGCM3")    
# GCMs <-  GCMs[1:2] #XXX for testing purposes

BGCcolors <- read.csv("C:\\Colin\\Projects\\Projects\\2019_CCISS\\InputData\\BGCzone_Colorscheme.csv")

# Knowledge Tables
treesuit="TreeSpp_ESuit_v11_18"
SiteSeries_Use <-read.csv(paste("InputData/","SiteSeries_Use_5",".csv",sep=""),stringsAsFactors=FALSE,na.strings=".")
spps.lookup <- read.csv("InputData\\Tree speciesand codes_2.0_2May2019.csv")

#===============================================================================
# create a dem from the climateBC input data
#===============================================================================

## create a dem from the climateBC input data
points <- read.csv(paste("InputData\\",grid,".csv", sep=""))
dim(points)

# points <- points[order(points$lon, points$lat),]
# head(points)

# transform latlon points to albers (this only works because the points were originally in albers)
P4S.latlon <- CRS("+proj=longlat +datum=WGS84")
P4S.albers <- CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs")
pts <- points[,4:3]
coordinates(pts) <- pts
projection(pts) <- P4S.latlon
pts.albers <- spTransform(pts, P4S.albers) #reproject the gridpoints

## assemble a simple spatial data frame of the points (seems awkward, but this is how i got the raster to work properly)
pts.df <- as.data.frame(pts.albers)
pts.df<- round(pts.df[,3:4],1)
df = data.frame(z = 1:dim(points)[1],
                xc = pts.df[,1],
                yc = pts.df[,2])
coordinates(df) = ~xc+yc
gridded(df) = TRUE
df = as(df, "SpatialGridDataFrame") # to full grid
plotOrder <- df$z #plotting order relative to input data

## rasterize
dem <- raster(df)
projection(dem) <- P4S.albers
par(mfrow=c(1,1))
plot(dem)
sum(!is.na(values(dem)))

X <- dem
values(X) <- points$el[df$z]
plot(X)

#===============================================================================
# generate the vector of mapped BGCs
#===============================================================================

BGC <- points$ID2
table(BGC)

BGC <- gsub(" ","",BGC)
table(BGC)
sort(table(BGC))

# merge small units
BGC[which(BGC=="CMAwh")] <- "CMAun"
BGC[which(BGC=="MHun")] <- "MHunp"
BGC[which(BGC=="ESSFxvw")] <- "ESSFxvp"
BGC[which(BGC=="ESSFdcp")] <- "ESSFdcw"

#BGC zones
BGCcolors <- read.csv("C:\\Colin\\Projects\\2019_CCISS\\InputData\\BGCzone_Colorscheme.csv")
zone <- rep(NA, length(BGC))
for(i in BGCcolors$zone){ zone[grep(i,BGC)] <- i }
table(zone)

#===============================================================================
# generic spatial data
#===============================================================================
### admin boundaries
bdy.bc <- readOGR("InputData\\BC_AB_US_Shp\\ProvincialOutline.shp")


#===============================================================================
# climate data for reference period
#===============================================================================

setwd("C:\\Colin\\Projects\\2019_CCISS")

fplot=paste("InputData\\", grid, "_Normal_1961_1990MSY.csv", sep="")

Columns = c("ID1", "ID2", "Latitude", "Longitude", "Elevation", "AHM", "bFFP",
            "CMD07","DD5_sp","EMT","Eref_sm","EXT","FFP","MCMT","MSP",
            "PPT07","PPT08", "PPT05","PPT06","PPT09", "SHM","TD","Tmax_sp","Tmin_at",
            "Tmin_sm","Tmin_wt", "PPT_at","PPT_wt", "PAS","eFFP",
            "Eref09","MAT","Tmin_sp","CMD")

Y0 <- fread(fplot, select = Columns, stringsAsFactors = FALSE, data.table = FALSE) #fread is faster than read.csv

Y0 <- Y0[!is.na(Y0[,2]),]

#####generate some additional variables
Y0$PPT_MJ <- Y0$PPT05 + Y0$PPT06 # MaY/June precip
Y0$PPT_JAS <- Y0$PPT07 + Y0$PPT08 + Y0$PPT09 # July/Aug/Sept precip
Y0$PPT.dormant <- Y0$PPT_at + Y0$PPT_wt # for calculating spring deficit
Y0$CMD.def <- 500 - (Y0$PPT.dormant)# start of growing season deficit original value was 400 but 500 seems better
Y0$CMD.def [Y0$CMD.def < 0] <- 0 #negative values set to zero = no deficit
Y0$CMDMax <- Y0$CMD07
Y0$CMD.total <- Y0$CMD.def + Y0$CMD

MAT.ref <- Y0$MAT
MSP.ref <- Y0$MSP
PPT_wt.ref <- Y0$PPT_wt
MAT.mean.ref <- mean(MAT.ref, na.rm=T)
MSP.mean.ref <- mean(MSP.ref, na.rm=T)
PPT_wt.mean.ref <- mean(PPT_wt.ref, na.rm=T)


#===============================================================================
# climate data for future periods
#===============================================================================
rcp="rcp45"
GCM=GCMs[1]
rcps <- c("rcp45", "rcp85")

  for(GCM in GCMs){
    Y1 <- fread(paste("InputData\\", grid, "_", GCM, "_BioVars.csv", sep=""), stringsAsFactors = FALSE, data.table = FALSE)
    Y1 <- Y1[!is.na(Y1[,2]),]
    
    #####generate some additional variables
    Y1$PPT_MJ <- Y1$PPT05 + Y1$PPT06 # MaY/June precip
    Y1$PPT_JAS <- Y1$PPT07 + Y1$PPT08 + Y1$PPT09 # July/Aug/Sept precip
    Y1$PPT.dormant <- Y1$PPT_at + Y1$PPT_wt # for calculating spring deficit
    Y1$CMD.def <- 500 - (Y1$PPT.dormant)# start of growing season deficit original value was 400 but 500 seems better
    Y1$CMD.def [Y1$CMD.def < 0] <- 0 #negative values set to zero = no deficit
    Y1$CMDMax <- Y1$CMD07
    Y1$CMD.total <- Y1$CMD.def + Y1$CMD
    
    ## assign single vectors to RCPs and proj.years
    Ystr <- strsplit(Y1[,1], "_")
    Y4 <- matrix(unlist(Ystr), ncol=3, byrow=TRUE)
    proj.years <- unique(Y4[,3])
    for(rcp in rcps){
      for(proj.year in proj.years){
      assign(paste("MAT", GCM, rcp, proj.year, sep="."), Y1$MAT[which(Y4[,2]==rcp & Y4[,3]==proj.year)])
      assign(paste("MSP", GCM, rcp, proj.year, sep="."), Y1$MSP[which(Y4[,2]==rcp & Y4[,3]==proj.year)])
      assign(paste("PPT_wt", GCM, rcp, proj.year, sep="."), Y1$PPT_wt[which(Y4[,2]==rcp & Y4[,3]==proj.year)])
    }
      print(rcp)
    }
    print(GCM)
  }


# calculate mean climate values for each GCM/year/rcp
for(rcp in rcps){
  for(proj.year in proj.years){
    MAT.change.mean <- rep(NA, length(GCMs))
    MSP.change.mean <- rep(NA, length(GCMs))
    PPT_wt.change.mean <- rep(NA, length(GCMs))
    for(GCM in GCMs){
      MAT <- get(paste("MAT", GCM, rcp, proj.year, sep="."))
      MSP <- get(paste("MSP", GCM, rcp, proj.year, sep="."))
      PPT_wt <- get(paste("PPT_wt", GCM, rcp, proj.year, sep="."))
      MAT.change <- MAT-MAT.mean.ref
      MSP.change <- MSP/MSP.mean.ref
      PPT_wt.change <- PPT_wt/PPT_wt.mean.ref
      MAT.change.mean[which(GCMs==GCM)] <- mean(MAT.change, na.rm=T)
      MSP.change.mean[which(GCMs==GCM)] <- mean(MSP.change, na.rm=T)
      PPT_wt.change.mean[which(GCMs==GCM)] <- mean(PPT_wt.change, na.rm=T)
    }
    assign(paste("MAT.change", rcp, proj.year, sep="."), MAT.change.mean)
    assign(paste("MSP.change", rcp, proj.year, sep="."), MSP.change.mean)
    assign(paste("PPT_wt.change", rcp, proj.year, sep="."), PPT_wt.change.mean)
    print(proj.year)
  }
  print(rcp)
}


#===============================================================================
# Spp suitability metrics relative to temperature change
#===============================================================================
edatopes<- c("B2", "C4", "D6")

edatope.name <- c("Subxeric-poor", "Mesic-medium", "Hygric-rich")
proj.year.name=c("2020s", "2050s", "2080s")
rcp.name=c("RCP4.5", "RCP8.5")
rcp=rcps[1]
proj.year=proj.years[2]
edatope="C4"

# Import suitability tables
wd="InputData"
treesuit2=paste(wd,"/",treesuit,".csv",sep="")
S1 <- read.csv(treesuit2,stringsAsFactors=F,na.strings=".")
S1 <- unique(S1)

## non-THLB BGCs for exclusion from results
BGCs_notin_THLB <- read.csv("InputData\\BGCs_notin_THLB.csv")
BGCs_notin_THLB <- BGCs_notin_THLB$BGC[which(BGCs_notin_THLB$Exlude=="x")]
exclude <- which(BGC%in%BGCs_notin_THLB)

# select the species to run the analysis on
spps <- unique(S1$Spp)
spps <- spps[-which(spps=="X")]
spps.candidate <- spps.lookup$TreeCode[which(spps.lookup$Exclude=="" & spps.lookup$Native=="N")]
spps.candidate <- spps.candidate[-which(spps.candidate=="Acb")]
spps <- spps[which(spps%in%spps.candidate)] 

# spps <- c("Ra", "Qg")
# spps <- c("Pl", "Fd", "Cw", "Sx", "At", "Py")
# spps <- c("Ba", "Bl", "Bg", "Yc", "Pa", "Hm", "Lw", "Hw","Dr", "Ep", "Act", "Sb", "Mb", "Ss")
# spps <- c("Pl", "Fd", "Cw", "Sx", "At", "Py", "Ba", "Bl", "Bg", "Yc", "Pa", 
#           "Hm", "Lw", "Hw","Dr", "Ep", "Act", "Sb", "Mb", "Ss", "Ra", "Qg")

# spps <- c("Pl", "Fd", "Cw")
# spps <- c("Ba", "Bl", "Bg")
# spps <- c("Yc", "Pa", "Hm")
# spps <- c("Lw", "Hw", "Py")
# spps <- c("Dr", "Ep", "At")

spp="Sx"
edatope=edatopes[2]
for(edatope in edatopes){
  for(spp in spps){
    
    Suit.ref <- read.csv(paste("OutputData\\Suit.ref", grid, spp, edatope, "csv", sep="."))[-exclude,1]
    # Suit.ref <- read.csv(paste("OutputData\\Suit.ref", grid, spp, edatope, "csv", sep="."))[,1]
    Suit.ref[Suit.ref==5] <- NA
    outRange.ref <- is.na(Suit.ref)
    Suit.ref[is.na(Suit.ref)] <- 5
    Suit.ref <- 1-(Suit.ref-1)/4
    
    for(hist.year in hist.years){
      Suit.proj <- read.csv(paste("OutputData\\Suit", grid, hist.year, spp, edatope, "csv", sep="."))[-exclude,1]
      # Suit.proj <- read.csv(paste("OutputData\\Suit", grid, hist.year, spp, edatope, "csv", sep="."))[,1]
      changeSuit <- Suit.ref-Suit.proj
      outRange <- outRange.ref
      outRange[which(changeSuit!=0)] <- FALSE
      # changeSuit[outRange==T] <- NA
      
      Suit.proj[is.na(Suit.proj)] <- 5
      Suit.proj <- 1-(Suit.proj-1)/4
      # hist(Suit.proj)
      
      assign(paste("SuitChange", hist.year, spp, edatope, sep="."), sum(Suit.proj)/sum(Suit.ref)) # note this is the sum of persistence and expansion
      assign(paste("persistence", hist.year, spp, edatope, sep="."), sum(Suit.proj[outRange.ref==F])/sum(Suit.ref))
      assign(paste("expansion", hist.year, spp, edatope, sep="."), sum(Suit.proj[outRange.ref==T])/sum(Suit.ref))
      assign(paste("expansionAbs", hist.year, spp, edatope, sep="."), sum(Suit.proj[outRange.ref==T]))
      # persistence+expansion
      # print(hist.year)
    }
    
    for(rcp in rcps){
      for(proj.year in proj.years){

        SuitChange <- rep(NA, length(GCMs))
        persistence <- rep(NA, length(GCMs))
        expansion <- rep(NA, length(GCMs))
        expansionAbs <- rep(NA, length(GCMs))

        for(GCM in GCMs){

          Suit.proj <- read.csv(paste("OutputData\\Suit", grid, GCM, rcp, proj.year, spp, edatope, "csv", sep="."))[-exclude,1]
          # Suit.proj <- read.csv(paste("OutputData\\Suit", grid, GCM, rcp, proj.year, spp, edatope, "csv", sep="."))[,1]
          changeSuit <- Suit.ref-Suit.proj
          outRange <- outRange.ref
          outRange[which(changeSuit!=0)] <- FALSE
          # changeSuit[outRange==T] <- NA

          Suit.proj[is.na(Suit.proj)] <- 5
          Suit.proj <- 1-(Suit.proj-1)/4
          # hist(Suit.proj)

          SuitChange[which(GCMs==GCM)] <- sum(Suit.proj)/sum(Suit.ref) # note this is the sum of persistence and expansion
          persistence[which(GCMs==GCM)] <- sum(Suit.proj[outRange.ref==F])/sum(Suit.ref)
          expansion[which(GCMs==GCM)] <- sum(Suit.proj[outRange.ref==T])/sum(Suit.ref)
          expansionAbs[which(GCMs==GCM)] <- sum(Suit.proj[outRange.ref==T])
          # persistence+expansion
          # print(GCM)
        }

        assign(paste("SuitChange", rcp, proj.year, spp, edatope, sep="."), SuitChange)
        assign(paste("persistence", rcp, proj.year, spp, edatope, sep="."), persistence)
        assign(paste("expansion", rcp, proj.year, spp, edatope, sep="."), expansion)
        assign(paste("expansionAbs", rcp, proj.year, spp, edatope, sep="."), expansionAbs)

        # print(proj.year)
      }
      # print(rcp)
    }
    
    print(spp)
  }
  print(edatope)
}


# Compile the metrics for all time periods/scenarios into a single vector

edatope=edatopes[2]
for(edatope in edatopes){
for(spp in spps){
  SuitChange <- vector()
  persistence <- vector()
  expansion <- vector()
  expansionAbs <- vector()
  for(rcp in rcps){
    for(proj.year in proj.years){
      SuitChange <- c(SuitChange, get(paste("SuitChange", rcp, proj.year, spp, edatope, sep=".")))
      persistence <- c(persistence, get(paste("persistence", rcp, proj.year, spp, edatope, sep=".")))
      expansion <- c(expansion, get(paste("expansion", rcp, proj.year, spp, edatope, sep=".")))
      expansionAbs <- c(expansionAbs, get(paste("expansionAbs", rcp, proj.year, spp, edatope, sep=".")))
      # print(proj.year)
    }
    # print(rcp)
  }

  assign(paste("SuitChange", spp, edatope, sep="."), SuitChange)
  assign(paste("persistence", spp, edatope, sep="."), persistence)
  assign(paste("expansion", spp, edatope, sep="."), expansion)
  assign(paste("expansionAbs", spp, edatope, sep="."), expansionAbs)
  
}
print(edatope)
}

# Compile the MAT for all time periods/scenarios into a single vector

MAT.change <- vector()
for(rcp in rcps){
  for(proj.year in proj.years){
    MAT.change <- c(MAT.change, get(paste("MAT.change", rcp, proj.year, sep=".")))
    # print(proj.year)
  }
  # print(rcp)
}


  
  ################
  ## manuscript figure
  ################
  library(msir)
  
    #vectors of projection specs. 
  seq.rcp <- NA
  seq.proj.year <-  NA
  seq.GCM <-  NA
  for(rcp in rcps){
    for(proj.year in proj.years){
      for(GCM in GCMs){
        seq.rcp <- c(seq.rcp, rcp)
        seq.proj.year <- c(seq.proj.year, proj.year)
        seq.GCM <- c(seq.GCM, GCM)
      }
    }
  }

  for(spp.focal in spps){
    
  # spp.focal <- "Hw"
  rcp="rcp45"
  proj.year=proj.years[2]
  # for(rcp in rcps){
  # for(proj.year in proj.years){
  png(filename=paste("Results\\Manu_BubblePlots\\CCISS.manu.ExpansionVsPersistence.THLB",spp.focal, proj.year, rcp,"png",sep="."), type="cairo", units="in", width=6.5, height=6.75, pointsize=13, res=400)
  mat <- matrix(c(1,2,3,4,4,4,5,5,6,5,5,7,9,9,9,8,8,8),6, byrow=T)   #define the plotting order
  layout(mat, widths=c(1,1,1), heights=c(0.8,0.15,1,1,0.07,0.1))   #set up the multipanel plot
  
  #===============================================================================
  # ABC. metric relative to temperature change
  
  ColScheme=c(2,1,3)
  
  metrics <- c("persistence", "expansion", "SuitChange")
  metric.names <- c("Persistence", "Expansion", "Suitable area (range)")
  par(mar=c(1.25,4,0.1,0.1), mgp=c(2.5,0.25,0))
  
  for(metric in metrics){
    numberDummies <- 15
    x <- c(rep(0,numberDummies),MAT.change)
    # ylim <- if(metric=="persistence") c(0,1.15) else c(-4,2)  #LOG SCALING OPTION
    ylim <- c(0,2.9)
    plot(0, xlim=c(0,max(MAT.change)+1.2), ylim=ylim, xaxs="i", col="white", xaxt="n", yaxt="n", 
         xlab="", ylab=metric.names[which(metrics==metric)], cex.lab=1.2)
    axis(1, at=0:7, labels = 0:7, tck=0)
    # if(metric=="persistence") {axis(2,at=seq(ylim[1], ylim[2], 0.2), labels=paste(seq(ylim[1], ylim[2], 0.2)*100,"%", sep=""), las=2, tck=0)   #LOG SCALING OPTION
    # } else axis(2,at=seq(ylim[1], ylim[2]), labels=paste(round(2^(seq(ylim[1], ylim[2]))*100,0),"%", sep=""), las=2, tck=0)
    axis(2,at=seq(ylim[1], ylim[2], 0.5), labels=paste(seq(ylim[1], ylim[2], 0.5)*100,"%", sep=""), las=2, tck=0)
    Common <- as.character(spps.lookup$EnglishName[which(spps.lookup$TreeCode==spp.focal)])
    Latin <- as.character(spps.lookup$ScientificName[which(spps.lookup$TreeCode==spp.focal)])
    panel <- paste("(", LETTERS[which(metrics==metric)],")", sep="")
    mtext(panel, side=3, line=-1.25, adj=0.05, cex=0.7, font=1)
    mtext(bquote(.(Common)~"("*.(spp.focal)*")"), side=3, line=-1.45, adj=0.955, cex=0.7, font=1)
    
    if(metric=="SuitChange") {
      x.focal <- x[numberDummies+which(seq.rcp==rcp & seq.proj.year==proj.year)]
      boxplot(x.focal, add=T, horizontal=TRUE, axes=FALSE, range=0, at=ylim[1]+0.025, boxwex = 0.12)
      text(max(x.focal), ylim[1]+0.025, paste(rcp.name[which(rcps==rcp)], ", ", proj.year.name[which(proj.years==proj.year)], sep=""), pos=4, cex=0.8)
    }
    
    # if(metric!="expansion") lines(c(-100, 99), if(metric=="persistence") c(1,1) else c(0,0), lty=2)
    # if(spp==spps.select[1])text(max(MAT.change)+1, 0.25, "1961-1990 baseline suitability", pos=2)
    
    for(edatope in edatopes[c(1,3,2)]){
      
      y <- c(if(metric=="expansion") rep(0,numberDummies) else rep(1,numberDummies),get(paste(metric, spp.focal, edatope, sep=".")))
      # if(metric!="persistence") {   #LOG SCALING OPTION
      #   y[y<2^(ylim[1])] <- 2^(ylim[1])
      #   y <- log2(y)
      # }
      # if(metric=="expansion") y[1] <- -4.5
      
      if(edatope=="C4") points(x,y, pch=16, col=ColScheme[which(edatopes==edatope)])
      l <- loess.sd(y~x, span=1, nsigma = 1)
      polygon(c(l$x, rev(l$x)), c(l$upper, rev(l$lower)), col=alpha(ColScheme[which(edatopes==edatope)], 0.5), border=NA)
      l <- loess(y~x)
      par(xpd=T)
      text(max(x)-0.25, predict(l, max(x)), edatope, pos=4, font= if(edatope=="C4") 2 else 1, col=ColScheme[which(edatopes==edatope)], cex= if(edatope=="C4") 1.1 else 1)
      par(xpd=F)
      
      print(edatope)
    }    
    print(metric)
  }
  
  par(mar=c(0,0,0,0))
  plot(1, type="n", axes=F, xlab="", ylab="")  
  text(1,1.25, bquote(Projected~temperature~change~relative~to~"1961-90"~"("*degree*C*")"), font=2,cex=1.2)  
  
  
  #===============================================================================
  # DEF. bubble plots
  
  colors = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)][-1]
  set.seed(2)
  ColScheme <- c(brewer.pal(n=12, "Paired"),sample(colors,length(spps)-12))
  
  
  for(edatope in edatopes[c(2,1,3)]){
    
    if(edatope=="C4") par(mar=c(0.1,3.5,0.1,0.1), mgp=c(2.5, 0.25, 0)) else par(mar=c(0.1,0.1,0.1,0.1), mgp=c(2.5, 0.25, 0))
      
    xlim <- c(0,1.05)
    ylim <- c(-5,3)
    plot(0, xlim=xlim, ylim=ylim, col="white", xaxt="n", yaxt="n", xlab="", ylab="")
    par(xpd=T)
    if(edatope!="B2") axis(1,at=seq(xlim[1], xlim[2], 0.2), labels=paste(seq(xlim[1], xlim[2], 0.2)*100,"%", sep=""), tck=0)
    par(xpd=F)
    if(edatope=="C4") axis(2,at=seq(ylim[1], ylim[2]), labels=paste(round(2^(seq(ylim[1], ylim[2]))*100),"%", sep=""), las=2, tck=0)
    title(ylab="Expansion beyond historically suitable range", cex.lab=1.2)
    iso <- seq(0,1.2, 0.001)
    lines(1-iso, log2(iso), lty=2, lwd=2, col="darkgray")
    if(edatope=="C4"){
      text(0.15, log2( 0.85), "Growing suit. range", pos=3, srt=-60/(ylim[2]/xlim[2])/2, col="darkgray", font=2)
      text(0.15, log2( 0.85), "Shrinking suit. range", pos=1, srt=-60/(ylim[2]/xlim[2])/2, col="darkgray", font=2)
    }
    panel <- paste("(", LETTERS[4:6][which(edatopes[c(2,1,3)]==edatope)],")", sep="")
    mtext(paste(panel," ", edatope," edatope", " (", edatope.name[which(edatopes==edatope)], " sites)", sep=""), side=3, line=-1.25, adj= if(edatope=="C4") 0.025 else 0.075, cex=0.7, font=1)
    
    
    for(spp in spps){
      
      x <- c(1,get(paste("persistence", spp, edatope, sep=".")))[which(seq.rcp==rcp & seq.proj.year==proj.year)]
      y <- c(0,get(paste("expansion", spp, edatope, sep=".")))[which(seq.rcp==rcp & seq.proj.year==proj.year)]
      y[y<2^(ylim[1])] <- 2^(ylim[1])
      y <- log2(y)
      
      # points(x,y)
      library(car)
      dataEllipse(x, y, levels=0.5, center.pch=21, add=T, col=ColScheme[which(spps==spp)], fill=T, lwd=0.5, plot.points=F)
      points(mean(x),mean(y), pch=21, bg=ColScheme[which(spps==spp)], cex=if(spp==spp.focal) 4.5 else 3)
      text(mean(x),mean(y), spp, cex=if(spp==spp.focal) 1 else 0.7, font=2)
      
      print(spp)
    }
    box()
  }
  
  par(mar=c(0,0,0,0))
  plot(1, type="n", axes=F, xlab="", ylab="")  
  text(1,1, "Persistence within historically suitable range",cex=1.2)  
  
  dev.off()
  
  print(proj.year)
  # }
  print(rcp)
  # }
}
  
  
  
  
  ################
  ## GIF
  ################

  library(tweenr)
  ## interpolate the time steps for an animated plot. 
  
  for(edatope in edatopes){
  for(spp in spps){
    rcp="rcp45"
    
  data <- data.frame(
    x = rnorm(length(GCMs), 1, 0.0001),
    y = rnorm(length(GCMs), 0, 0.0001),
    z = rnorm(length(GCMs), 0, 0.0001),
    time = 1,
    group = GCMs,
    ease = rep('sine-in-out', length(GCMs))
  )
  
  # for(hist.year in c(1995, 2005, 2014)){
    for(hist.year in c(2004)){
      data <- rbind(data, data.frame(
      x = rnorm(length(GCMs), get(paste("persistence", hist.year, spp, edatope, sep=".")), 0.0001),
      y = rnorm(length(GCMs), get(paste("expansion", hist.year, spp, edatope, sep=".")), 0.0001),
      z = rnorm(length(GCMs), 0, 0.0001),
      time = max(data$time)+1,
      group = GCMs,
      ease = rep('sine-in-out', length(GCMs))
    ))
    
    # print(hist.year)
    
  }
  
  for(proj.year in c(proj.years)){

    data <- rbind(data, data.frame(
      x = get(paste("persistence", rcp, proj.year, spp, edatope, sep=".")),
      y = get(paste("expansion", rcp, proj.year, spp, edatope, sep=".")),
      z = rep(0, length(GCMs)),
      time = max(data$time)+1,
      group = GCMs,
      ease = rep('sine-in-out', length(GCMs))
    ))
    
    # print(proj.year)
  }
  
  rcp="rcp85"
    data <- rbind(data, data.frame(
    x = get(paste("persistence", rcp, proj.year, spp, edatope, sep=".")),
    y = get(paste("expansion", rcp, proj.year, spp, edatope, sep=".")),
    z = rep(0, length(GCMs)),
    time = max(data$time)+1,
    group = GCMs,
    ease = rep('sine-in-out', length(GCMs))
  ))
  
  nframes=(max(data$time)-1)*50+1
  assign(paste("data", spp, edatope, sep="."), tween_elements(data, 'time', 'group', 'ease', nframes = nframes))
  
  print(spp)
  }
  print(edatope)
  }
  
  ######
  ## plot the individual frames
  ######

  edatope="C4"
  
  
    colors = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)][-1]
  set.seed(2)
  ColScheme <- c(brewer.pal(n=12, "Paired"),sample(colors,length(spps)-12))
  
  spp.focal="x"
  
  
  # for(edatope in edatopes[c(2,1,3)]){
  for(j in 1:nframes){
  
  png(filename=paste("Results\\BubbleGIF\\ExpansionVsPersistence",100+j,"png",sep="."), type="cairo", units="in", width=6.5, height=6.5, pointsize=13, res=300)
  par(mar=c(3,4.5,0.1,0.1), mgp=c(2.75, 0.25, 0), xpd=F)
  
  xlim <- c(0,1.05)
  ylim <- c(-5,3)
  plot(0, xlim=xlim, ylim=ylim, col="white", xaxt="n", yaxt="n", xlab="", ylab="")
  axis(1,at=seq(xlim[1], xlim[2], 0.2), labels=paste(seq(xlim[1], xlim[2], 0.2)*100,"%", sep=""), tck=0)
  axis(2,at=seq(ylim[1], ylim[2]), labels=paste(round(2^(seq(ylim[1], ylim[2]))*100),"%", sep=""), las=2, tck=0)
  title(ylab="Expansion beyond historically suitable range", cex.lab=1.2)
  par(mgp=c(1.5, 0.25, 0))
  title(xlab="Persistence within historically suitable range", cex.lab=1.2)
  iso <- seq(0,1.2, 0.001)
  lines(1-iso, log2(iso), lty=2, lwd=2, col="darkgray")
  if(edatope=="C4"){
    text(0.15, log2( 0.85), "Growing suit. range", pos=3, srt=-60/(ylim[2]/xlim[2])/2, col="darkgray", font=2)
    text(0.15, log2( 0.85), "Shrinking suit. range", pos=1, srt=-60/(ylim[2]/xlim[2])/2, col="darkgray", font=2)
  }
  panel <- paste("(", LETTERS[4:6][which(edatopes[c(2,1,3)]==edatope)],")", sep="")
  # mtext(paste(panel," ", edatope," edatope", " (", edatope.name[which(edatopes==edatope)], " sites)", sep=""), side=3, line=-1.25, adj= if(edatope=="C4") 0.025 else 0.075, cex=0.7, font=1)
  
  for(spp in spps){
    
    data <- get(paste("data", spp, edatope, sep="."))
    x <- data$x[which(data$.frame==j)]
    y <- data$y[which(data$.frame==j)]
    y[y<2^(ylim[1])] <- 2^(ylim[1])
    y <- log2(y)
    
    # points(x,y)
    library(car)
    dataEllipse(x, y, levels=0.5, center.pch=21, add=T, col=ColScheme[which(spps==spp)], fill=T, lwd=0.5, plot.points=F)
    points(mean(x),mean(y), pch=21, bg=ColScheme[which(spps==spp)], cex=if(spp==spp.focal) 4.5 else 3)
    text(mean(x),mean(y), spp, cex=if(spp==spp.focal) 1 else 0.7, font=2)
    
    # print(spp)
  }
  box()
  if(j%in%c(0:20, 40:60,90:110,140:160,190:210,240:251)) mtext(rep(c("1970s", "2000s", "2020s (RCP4.5)", "2050s (RCP4.5)", "2080s (RCP4.5)", "2080s (RCP8.5)"),each=50)[j+25], 3, -1.5, adj=0.975, font=2)
  # print(edatope)
  #   }
  
  dev.off()
  
  print(j)
  }
  
  # setwd("C:\\Colin\\Projects\\2019_CCISS\\Results\\BubbleGIF")
  # ## compile the gif
  # library(magick)
  # library(animation)
  # ani.options(convert = 'C:\\Program Files\\ImageMagick-7.0.8-Q16\\magick.exe')
  # system("magick -loop 1 -delay 0 *.png test.gif")
  
  setwd("C:\\Colin\\Projects\\2019_CCISS\\Results\\BubbleGIF_2000s")
  ## compile the gif
  library(magick)
  ani.options(convert = 'C:\\Program Files\\ImageMagick-7.0.8-Q16\\magick.exe')
  library(animation)
  system("magick -loop 1 -delay 0 *.png test.gif")
  
  setwd("C:\\Colin\\Projects\\2019_CCISS\\Results\\BubbleGIF_2020s")
  ## compile the gif
  library(magick)
  ani.options(convert = 'C:\\Program Files\\ImageMagick-7.0.8-Q16\\magick.exe')
  library(animation)
  system("magick -loop 1 -delay 0 *.png test.gif")
  
  setwd("C:\\Colin\\Projects\\2019_CCISS\\Results\\BubbleGIF_2050s")
  ## compile the gif
  library(magick)
  ani.options(convert = 'C:\\Program Files\\ImageMagick-7.0.8-Q16\\magick.exe')
  library(animation)
  system("magick -loop 1 -delay 0 *.png test.gif")
  
  setwd("C:\\Colin\\Projects\\2019_CCISS\\Results\\BubbleGIF_2080s")
  ## compile the gif
  library(magick)
  library(animation)
  ani.options(convert = 'C:\\Program Files\\ImageMagick-7.0.8-Q16\\magick.exe')
  system("magick -loop 1 -delay 0 *.png test.gif")
  
  setwd("C:\\Colin\\Projects\\2019_CCISS\\Results\\BubbleGIF_RCP85")
  ## compile the gif
  library(magick)
  ani.options(convert = 'C:\\Program Files\\ImageMagick-7.0.8-Q16\\magick.exe')
  library(animation)
  system("magick -loop 1 -delay 0 *.png test.gif")
  
  setwd("C:\\Colin\\Projects\\2019_CCISS")
  

  
  
  
  
  
  