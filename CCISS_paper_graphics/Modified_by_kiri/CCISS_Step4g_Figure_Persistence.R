
##======================================================================================
## CCISS Publication Scripts
## Step 4g - Figures - Manuscript plot of Suitability and species persistence
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

GCMs <-  c("ACCESS1-0","CanESM2","CCSM4","CESM1-CAM5","CNRM-CM5","CSIRO-Mk3-6-0", "GFDL-CM3","GISS-E2R", "HadGEM2-ES", "INM-CM4", "IPSL-CM5A-MR", "MIROC-ESM", "MIROC5", "MPI-ESM-LR","MRI-CGCM3")
rcps <- c("rcp45", "rcp85")
proj.years <- c(2025, 2055, 2085)
hist.years <- c(1995, 2004, 2005, 2009, 2014, 2017)
edatopes<- c("B2", "C4", "D6")
edatope.name <- c("nutrient poor, subxeric", "nutrient-medium, mesic", "nutrient-rich, hygric")
proj.year.name=c("2020s", "2050s", "2080s")
rcp.name=c("RCP4.5", "RCP8.5")

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
# generate the vector of predicted referencc BGCs 
#===============================================================================

BGC <- as.character(read.csv(paste("OutputData\\BGC.pred", grid, "ref.csv", sep="."))[,1])

#BGC zones
BGCcolors <- read.csv("C:\\Colin\\Projects\\2019_CCISS\\InputData\\BGCzone_Colorscheme.csv")
zone <- rep(NA, length(BGC))
for(i in BGCcolors$zone){ zone[grep(i,BGC)] <- i }
table(zone)
zone <- factor(zone, levels=c("CDF", "CWH", "MH", "ESSF", "MS", "IDF", "PP", "BG", "ICH", "SBPS", "SBS", "BWBS", "SWB", "CMA", "IMA", "BAFA"))


#===============================================================================
# generic spatial data
#===============================================================================
### admin boundaries
bdy.bc <- readOGR("InputData\\BC_AB_US_Shp\\ProvincialOutline.shp")


#===============================================================================
# Import suitability tables
#===============================================================================
wd="InputData"
treesuit2=paste(wd,"/",treesuit,".csv",sep="")
S1 <- read.csv(treesuit2,stringsAsFactors=F,na.strings=".")
S1 <- unique(S1)

S1[grep("CRF", S1$BGC),]
S1[which(S1$BGC=="CRFdhz"),]
S1[which(S1$BGC=="CWHxmz"),]

# select the species to run the analysis on
spps <- unique(S1$Spp)
spps <- spps[-which(spps=="X")]
spps.candidate <- spps.lookup$TreeCode[-which(spps.lookup$Exclude=="x")]
spps <- spps[which(spps%in%spps.candidate)] 

#===============================================================================
# calculate mean optionality and turnover for each edatope/year/rcp
#===============================================================================

# read in metrics for reference period and historical decades
for(edatope in edatopes){
  assign(paste("SuitRichness.ref", edatope, sep="."), read.csv(paste("OutputData\\SuitRichness.ref", grid, edatope, "csv", sep="."))[,1])
  assign(paste("SppRichness.ref", edatope, sep="."), read.csv(paste("OutputData\\SppRichness.ref", grid, edatope, "csv", sep="."))[,1])
  print(edatope)
}

for(edatope in edatopes){
  for(rcp in rcps){
    for(proj.year in proj.years){
      assign(paste("SuitRichness", rcp, proj.year, edatope, sep="."), read.csv(paste("OutputData\\SuitRichness", grid, rcp, proj.year, edatope, "csv", sep="."))[,1])
      assign(paste("SuitRichnessChange", rcp, proj.year, edatope, sep="."), read.csv(paste("OutputData\\SuitRichnessChange", grid, rcp, proj.year, edatope, "csv", sep="."))[,1])
      assign(paste("SppRichness", rcp, proj.year, edatope, sep="."), read.csv(paste("OutputData\\SppRichness", grid, rcp, proj.year, edatope, "csv", sep="."))[,1])
      assign(paste("SppRichnessChange", rcp, proj.year, edatope, sep="."), read.csv(paste("OutputData\\SppRichnessChange", grid, rcp, proj.year, edatope, "csv", sep="."))[,1])
      assign(paste("SuitTurnover", rcp, proj.year, edatope, sep="."), read.csv(paste("OutputData\\SuitTurnover", grid, rcp, proj.year, edatope, "csv", sep="."))[,1])
      assign(paste("SppTurnover", rcp, proj.year, edatope, sep="."), read.csv(paste("OutputData\\SppTurnover", grid, rcp, proj.year, edatope, "csv", sep="."))[,1])
      assign(paste("SuitPersistence", rcp, proj.year, edatope, sep="."), read.csv(paste("OutputData\\SuitPersistence", grid, rcp, proj.year, edatope, "csv", sep="."))[,1])
      assign(paste("SppPersistence", rcp, proj.year, edatope, sep="."), read.csv(paste("OutputData\\SppPersistence", grid, rcp, proj.year, edatope, "csv", sep="."))[,1])
      print(proj.year)
    }
    print(rcp)
  }
  print(edatope)
}



#===============================================================================
# BGC projections
#===============================================================================

# Future BGC
PredSum <- data.frame()
for(rcp in rcps){
  for(proj.year in proj.years){
    for(GCM in GCMs){
      BGC.pred <- as.character(read.csv(paste("OutputData\\BGC.pred",grid, GCM, rcp, proj.year,".csv", sep=""))[,1])
      assign(paste("BGC.pred", GCM, rcp, proj.year, sep="."), BGC.pred) #bgc projection
      PredSum <- rbind(PredSum, as.data.frame(table(BGC.pred)))
      # print(GCM)
    }
    print(proj.year)
  }
  print(rcp)
}


#===============================================================================
# calculate mean MAT change for each model prediction
#===============================================================================

setwd("C:\\Colin\\Projects\\2019_CCISS")

fplot=paste("InputData\\", grid, "_Normal_1961_1990MSY.csv", sep="")
Y0 <- fread(fplot, select = "MAT", stringsAsFactors = FALSE, data.table = FALSE) #fread is faster than read.csv
MAT.ref <- Y0$MAT
MAT.mean.ref <- mean(MAT.ref, na.rm=T)

for(hist.year in hist.years){
  Y0 <- fread(paste("InputData\\", grid, "_", hist.year, "_BioVars.csv", sep=""), select = "MAT", stringsAsFactors = FALSE, data.table = FALSE) #fread is faster than read.csv
  assign(paste("MAT", hist.year, sep="."), Y0$MAT)
  assign(paste("MAT.change", hist.year, sep="."), mean(Y0$MAT, na.rm=T)-MAT.mean.ref)
  print(hist.year)
}

for(GCM in GCMs){
  Y1 <- fread(paste("InputData\\", grid, "_", GCM, "_BioVars.csv", sep=""), select = c("GCM", "MAT"), stringsAsFactors = FALSE, data.table = FALSE)
  ## assign single vectors to RCPs and proj.years
  Ystr <- strsplit(Y1[,1], "_")
  Y4 <- matrix(unlist(Ystr), ncol=3, byrow=TRUE)
  for(rcp in rcps){
    for(proj.year in proj.years){
      assign(paste("MAT", GCM, rcp, proj.year, sep="."), Y1$MAT[which(Y4[,2]==rcp & Y4[,3]==proj.year)])
    }
    print(rcp)
  }
  print(GCM)
}

# calculate mean climate values for each GCM/year/rcp
for(rcp in rcps){
  for(proj.year in proj.years){
    MAT.change.mean <- rep(NA, length(GCMs))
    for(GCM in GCMs){
      MAT <- get(paste("MAT", GCM, rcp, proj.year, sep="."))
      MAT.change <- MAT-MAT.mean.ref
      MAT.change.mean[which(GCMs==GCM)] <- mean(MAT.change, na.rm=T)
    }
    assign(paste("MAT.change", rcp, proj.year, sep="."), MAT.change.mean)
    print(proj.year)
  }
  print(rcp)
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


#===============================================================================
# calculate mean persistence for each GCM/year/rcp
#===============================================================================

for(edatope in edatopes){
  for(hist.year in hist.years){
    SuitPersistence <- read.csv(paste("OutputData\\SuitPersistence", grid, hist.year, edatope, "csv", sep="."))[,1]
    SuitPersistence[SuitPersistence==99] <- NA # remove grid cells with no current suitability
    SuitPersistence.mean <- mean(SuitPersistence, na.rm=T)
    SppPersistence <- read.csv(paste("OutputData\\SppPersistence", grid, hist.year, edatope, "csv", sep="."))[,1]
    SppPersistence[SppPersistence==99] <- NA # remove grid cells with no current suitability
    SppPersistence.mean <- mean(SppPersistence, na.rm=T)
    assign(paste("SuitPersistence.mean", hist.year, edatope, sep="."), SuitPersistence.mean)
    assign(paste("SppPersistence.mean", hist.year, edatope, sep="."), SppPersistence.mean)
  }
  print(edatope)
}


for(edatope in edatopes){
  for(rcp in rcps){
    for(proj.year in proj.years){
      SuitPersistence.mean <- rep(NA, length(GCMs))
      SppPersistence.mean <- rep(NA, length(GCMs))
      for(GCM in GCMs){
        SuitPersistence <- read.csv(paste("OutputData\\SuitPersistence", grid, GCM, rcp, proj.year, edatope, "csv", sep="."))[,1]
        SuitPersistence[SuitPersistence==99] <- NA # remove grid cells with no current suitability
        SuitPersistence.mean[which(GCMs==GCM)] <- mean(SuitPersistence, na.rm=T)
        SppPersistence <- read.csv(paste("OutputData\\SppPersistence", grid, GCM, rcp, proj.year, edatope, "csv", sep="."))[,1]
        SppPersistence[SppPersistence==99] <- NA # remove grid cells with no current suitability
        SppPersistence.mean[which(GCMs==GCM)] <- mean(SppPersistence, na.rm=T)
      }
      assign(paste("SuitPersistence.mean", rcp, proj.year, edatope, sep="."), SuitPersistence.mean)
      assign(paste("SppPersistence.mean", rcp, proj.year, edatope, sep="."), SppPersistence.mean)
      print(proj.year)
    }
    print(rcp)
  }
  print(edatope)
}


# Compile the metrics for all time periods/scenarios into a single vector
edatope=edatopes[2]
for(edatope in edatopes){
  for(spp in spps){
    SuitPersistence.mean.all <- vector()
    SppPersistence.mean.all <- vector()
    for(rcp in rcps){
      for(proj.year in proj.years){
        SuitPersistence.mean.all <- c(SuitPersistence.mean.all, get(paste("SuitPersistence.mean", rcp, proj.year, edatope, sep=".")))
        SppPersistence.mean.all <- c(SppPersistence.mean.all, get(paste("SppPersistence.mean", rcp, proj.year, edatope, sep=".")))
        # print(proj.year)
      }
      # print(rcp)
    }
    assign(paste("SuitPersistence.mean.all", edatope, sep="."), SuitPersistence.mean.all)
    assign(paste("SppPersistence.mean.all", edatope, sep="."), SppPersistence.mean.all)
  }
  print(edatope)
}

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




################################
## Manuscript figure
###################################

edatope=edatopes[2]
rcp=rcps[1]
proj.year=proj.years[2]

for(proj.year in proj.years){
  for(edatope in edatopes){
    SuitRichness.ref <- get(paste("SuitRichness.ref", edatope, sep="."))
    SuitRichness.proj <- get(paste("SuitRichness", rcp, proj.year, edatope, sep="."))
    SuitRichnessChange <- get(paste("SuitRichnessChange", rcp, proj.year, edatope, sep="."))
    SuitRichnessChangePct <- SuitRichness.proj/SuitRichness.ref
    SuitRichnessChangePct[!is.finite(SuitRichnessChangePct)] <- NA
    SppRichness.ref <- get(paste("SppRichness.ref", edatope, sep="."))
    SppRichness.proj <- get(paste("SppRichness", rcp, proj.year, edatope, sep="."))
    SppRichnessChange <- get(paste("SppRichnessChange", rcp, proj.year, edatope, sep="."))
    SppRichnessChangePct <- SppRichness.proj/SppRichness.ref
    SppRichnessChangePct[!is.finite(SppRichnessChangePct)] <- NA
    SuitTurnover <- get(paste("SuitTurnover", rcp, proj.year, edatope, sep="."))
    SppTurnover <- get(paste("SppTurnover", rcp, proj.year, edatope, sep="."))
    SuitTurnover[which(SuitTurnover==99)] <- NA
    SppTurnover[which(SppTurnover==99)] <- NA
    SuitPersistence <- get(paste("SuitPersistence", rcp, proj.year, edatope, sep="."))
    SppPersistence <- get(paste("SppPersistence", rcp, proj.year, edatope, sep="."))
    
    metric <- "SppPersistence"
    
    # x11(width=6.5, height=5, pointsize=8)
    png(filename=paste("Results\\Manu_Persistence\\CCISS_manu_", metric, edatope, rcp, proj.year,"png",sep="."), type="cairo", units="in", width=6.5, height=4.55, pointsize=8, res=400)
    # pdf(file=paste("Results\\CCISS_SummaryByBGC_", metric,".pdf",sep=""),  width=7.5, height=5.625, pointsize=15)
    
    ylim <- c(0,1.4)
    y <- get(metric)
    # y[y<2^(ylim[1])] <- 2^(ylim[1])
    # y <- log2(y)
    
    par(mar=c(0,0,0,0))
    plot(0, xlim=c(0,1), ylim=c(0,1), col="white", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
    # box()
    
    #map
    par(mar=c(0,0,0,0), plt = c(0.05, 0.85, 0.025, 0.999), new = TRUE)
    # SuitTurnover <- get(paste("SuitTurnover", rcp, proj.year, edatope, sep="."))
    # SppTurnover <- get(paste("SppTurnover", rcp, proj.year, edatope, sep="."))
    values(X) <- y[plotOrder]
    
    breakpoints <- seq(0,2,0.2); length(breakpoints)
    labels <- c("0%", "50%", "No change", "150%", "200%")
    ColScheme <- c(brewer.pal(11,"RdBu")[1:4], rep("gray90", 2), brewer.pal(11,"RdBu")[8:11]); length(ColScheme)
    # ColScheme <- brewer.pal(11,"RdBu")[-6]; length(ColScheme)
    
    plot(bdy.bc, border="black", lwd=0.4)
    image(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks=breakpoints, maxpixels= ncell(X))
    mtext("(A)", side=3, line=-2.5, adj=0.0, cex=1, font=2)
    par(xpd=T)
    xl <- 250000; yb <- 1000000; xr <- 350000; yt <- 1500000
    # xl <- 2025000; yb <- 500000; xr <- 2100000; yt <- 950000
    rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme)
    text(rep(xr,length(labels)),seq(yb,yt,(yt-yb)/(length(labels)-1)),labels,pos=4,cex=1,font=1)
    # text(rep(xl,2),c(yb,yt)+c(-10000, 10000),c("No change", "Full turnover"),pos=c(1,3),cex=1,font=1)
    text(xl-40000, mean(c(yb,yt))-30000, paste("Suitability persistence\n(", proj.year.name[which(proj.years==proj.year)], ", ", c("RCP4.5", "RCP8.5")[which(rcps==rcp)], ")", sep=""), srt=90, pos=3, cex=1, font=2)
    # rect(xl,  yt+20000,  xr,  yt+60000,  col=ColScheme[length(ColScheme)])
    # text(xr,  yt+40000,  bquote(">"*.(breakseq[3])*sigma),pos=4,cex=1,font=1)  
    par(xpd=F)
    mtext(paste(edatope, " edatope (", edatope.name[which(edatopes==edatope)], " sites)", sep=""), side=3, line=-1.5, adj=0.35, cex=1, font=2)
    
    ##==============================
    ## Callout box for Location 1
    bgc.select <- "MSdk"
    q.select <- 0.9 # 0.90 is good
    focal.bgc <- which(BGC==bgc.select)
    # pt <- focal.bgc[which(SuitTurnover[focal.bgc]==min(SuitTurnover[focal.bgc][SuitTurnover[focal.bgc]>quantile(SuitTurnover[focal.bgc], q.select, na.rm=T)], na.rm=T))[1]] #select the grid point
    # pt <- pt[1]
    pt <- 96923
    lines(c(pts.df[pt,1], pts.df[pt,1]+200000), c(pts.df[pt,2], pts.df[pt,2]+0), lwd=1.5)
    points(pts.df[pt,], pch=21, cex=1.75, lwd=1.5, bg=ColScheme[cut(y[pt],breaks=breakpoints)])
    
    # assemble the species mix for the reference period
    comm.ref.pt <- rep(NA, length(spps))
    for(spp in spps){
      Suit <- read.csv(paste("OutputData\\Suit.ref", grid, spp, edatope, "csv", sep="."))[pt,1]
      Suit[is.na(Suit)] <- 5  #XXX note this is different from the equivalent line for the other time periods.
      Suit <- 1-(Suit-1)/4
      comm.ref.pt[which(spps==spp)] <- Suit
    }
    names(comm.ref.pt) <- spps
    #convert fractional suitability to standard suitability
    comm.ref.pt <- round((1-comm.ref.pt)*4+1)
    comm.ref.pt[comm.ref.pt==4] <- 5
    
    # assemble the projected community
    for(GCM in GCMs){
      comm.proj <- rep(NA, length(spps))
      for(spp in spps){
        Suit <- read.csv(paste("OutputData\\Suit", grid, GCM, rcp, proj.year, spp, edatope, "csv", sep="."))[pt,1]
        Suit[Suit==5] <- 5
        Suit <- 1-(Suit-1)/4
        comm.proj[which(spps==spp)] <- Suit
      }
      names(comm.proj) <- spps
      if(GCM==GCMs[1]) comm.proj.pt <- comm.proj else comm.proj.pt <- rbind(comm.proj.pt,comm.proj)
      # print(GCM)
    }
    comm.proj.pt.ensMean <- apply(comm.proj.pt, 2, FUN=mean, na.rm=T)
    comm.proj.pt.ensMean <- round((1-comm.proj.pt.ensMean)*4+1)
    comm.proj.pt.ensMean[comm.proj.pt.ensMean==4] <- 5
    
    std.ref <- paste(paste(names(comm.ref.pt)[which(comm.ref.pt==1)],collapse=""), "(", 
                     paste(names(comm.ref.pt)[which(comm.ref.pt==2)],collapse=""), ")((", 
                     paste(names(comm.ref.pt)[which(comm.ref.pt==3)],collapse=""), "))", sep="")
    std.proj <- paste(paste(names(comm.proj.pt.ensMean)[which(comm.proj.pt.ensMean==1)],collapse=""), "(", 
                      paste(names(comm.proj.pt.ensMean)[which(comm.proj.pt.ensMean==2)],collapse=""), ")((", 
                      paste(names(comm.proj.pt.ensMean)[which(comm.proj.pt.ensMean==3)],collapse=""), "))", sep="")
    
    BGC.list <- vector()
    for(GCM in GCMs){
      BGC.list[which(GCMs==GCM)] <- get(paste("BGC.pred", GCM, rcp, proj.year, sep="."))[pt]
      # print(BGC.pred[pt])
      # print(GCM)
    }
    BGC.list 
    
    par(mar=c(0,0,0,0), plt = c(0.825, 0.995, 0.01, 0.445), new = TRUE, mgp=c(2,0.1,0))
    plot(0, xlim=c(0,1), ylim=c(0,1), col="white", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
    rect(-9,-9,9,9, col="white", lwd=0)
    mtext("(D)", side=3, line=-1, adj=-0.2, cex=1, font=2)
    box()
    
    # par(mar=c(0,0,0,0), plt = c(0.6, 0.99, 0.01, 0.175), new = TRUE, mgp=c(2,0.1,0))
    # plot(0, xlim=c(0,1), ylim=c(0,1), col="white", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
    pos1 <- 1
    text(-0.025, pos1-0, paste("Historical bioclimate"), pos=4, cex=0.8, font=2, offset=0.1)
    text(0, pos1-0.05, bgc.select, pos=4, cex=0.8, offset=0.1)
    text(-0.025, pos1-.125, paste("Projected bioclimates"), pos=4, cex=0.8, font=2, offset=0.1)
    pos2 <- 0.44
    text(-.025, pos2-0, "Historical suitabilities", pos=4, cex=0.8, font=2, offset=0.1)
    text(0.0, pos2-0.05, std.ref, pos=4, cex=0.8, offset=0.1)
    text(-.025, pos2-0.125, paste("Projected suitabilities", sep=""), pos=4, cex=0.8, font=2, offset=0.1)
    for(i in 1:length(unique(BGC.list))){
      BGC.i <- names(rev(sort(table(BGC.list))))[i]
      comm.proj.pt.i <- comm.proj.pt[which(BGC.list==BGC.i)[1],]
      std.proj.i <- paste(paste(names(comm.proj.pt.i)[which(comm.proj.pt.i==1)],collapse=""), "(", 
                          paste(names(comm.proj.pt.i)[which(comm.proj.pt.i==0.75)],collapse=""), ")((", 
                          paste(names(comm.proj.pt.i)[which(comm.proj.pt.i==0.5)],collapse=""), "))", sep="")
      text(0.0, pos2-0.125-0.05*i, paste(BGC.i, ": ", std.proj.i, sep=""), pos=4, cex=0.75, offset=0.1)
    }
    
    text(-.025, 0, paste("Mean persistence: ", round(SuitPersistence[pt]*100), "%", sep=""), pos=4, cex=0.8, offset=0.1, font=2)
    
    par(mar=c(3,0,1,0), plt = c(0.855, 0.99, 0.29, 0.355), new = TRUE, mgp=c(0.5,0.2,0))
    bar.bgc <- names(rev(sort(table(BGC.list))))
    bar.zone <- rep(NA, length(bar.bgc))
    for(i in BGCcolors$zone){ bar.zone[grep(i,bar.bgc)] <- i }
    barplot(rev(sort(table(BGC.list))), col=as.character(BGCcolors$HEX[match(bar.zone, BGCcolors$zone)]), horiz=F, las=2, ylab=list("# GCMs", cex=0.7), yaxt="n", cex.names=0.7)
    par(mgp=c(1,0,0))
    axis(2, at=seq(0,12,4), labels=seq(0,12,4), las=1, tck=0, cex.axis=0.8, lwd=0)
    
    # spp.ref <- c(names(comm.ref.pt)[which(comm.ref.pt==1)], names(comm.ref.pt)[which(comm.ref.pt==2)], names(comm.ref.pt)[which(comm.ref.pt==3)])
    # boxplot(0, col="white", xlim=c(1,length(spp.ref)), ylim=c(0,1))
    # for(spp in spp.ref){
    #   y <- comm.proj.pt[,which(names(comm.proj.pt)==spp)]
    #   boxplot(y, add=T, at=which(spp.ref==spp), range=0, names=spp)
    # }
    # axis(1, at=1:length(spp.ref), labels = spp.ref, tck=0)
    
    
    ##===========================
    ## Summary by zone
    par(mar=c(4.5,2,0.1,0.1), plt = c(0.1, 0.325, 0.125, 0.325), new = TRUE, mgp=c(2.5,0.25,0))
    xlim=c(1, length(levels(droplevels(zone[-grep("BAFA|CMA|IMA",zone)]))))
    z <- boxplot(y~zone, ylim=ylim, ylab="", vertical = TRUE, plot=F)
    for(i in 1:length(levels(zone))){ 
      temp <- y[which(zone==levels(zone)[i])]
      z$stats[c(1,5), i] <- quantile(temp[!is.na(temp)],c(0.05, 0.95))
    }
    bxp(z, ylim=ylim, xlim=xlim, xaxt="n", yaxt="n", xaxs="i", ylab="", pch=0,outline=FALSE)
    lines(c(-99,99), c(1,1), lwd=2, col="darkgrey")
    bxp(z, add=T, boxfill = as.character(BGCcolors$HEX[match(levels(droplevels(zone[-grep("BAFA|CMA|IMA",zone)])), BGCcolors$zone)]), xaxt="n", yaxt="n",xaxs="i", 
        ylab=paste("Suitability persistence\n(", proj.year.name[which(proj.years==proj.year)], ", ", c("RCP4.5", "RCP8.5")[which(rcps==rcp)], ")", sep=""),  pch=0,outline=FALSE)
    axis(1, at=1:length(levels(zone)), levels(zone), tick=F, las=2, cex=0.8)
    axis(2, seq(ylim[1],ylim[2],.25), paste(seq(ylim[1],ylim[2],.25)*100, "%", sep=""), tick=F, las=2)
    # lines(range(grep("CWH|MH|CDF", levels(zone)))+c(-0.6,0.6), if(scenario=="6variable") rep(2.7,2) else rep(3.2,2), col=alpha("red", 0.5), lwd=3, lty=1)
    # lines(range(grep("ESSF|MS|IDF|PP|BG|ICH|SBPS|SBS|BWBS|SWB", levels(zone)))+c(-0.6,0.6), if(scenario=="6variable") rep(1.5,2) else rep(1.6,2), col=alpha("red", 0.5), lwd=3, lty=1)
    axis(1, at=c(mean(grep("CWH|MH|CDF", levels(zone))), mean(grep("ESSF|MS|IDF|PP|BG|ICH|SBPS|SBS|BWBS|SWB", levels(zone))),mean(grep("CMA|IMA|BAFA", levels(zone)))), 
         labels = c("Coast", "Interior", "Alpine"), tick=F, line=2.95 )
    par(xpd=NA)
    x1 <- range(grep("CWH|MH|CDF", levels(zone)))+c(-0.4,0.4)
    x2 <- range(grep("ESSF|MS|IDF|PP|BG|ICH|SBPS|SBS|BWBS|SWB", levels(zone)))+c(-0.4,0.4)
    # x3 <- range(grep("CMA|IMA|BAFA", levels(zone)))+c(-0.4,0.4)
    bracketpos=-0.7
    lines(rep(x1, each=2), bracketpos*c(0.95, 1,1,0.95))
    lines(rep(x2, each=2), bracketpos*c(0.95, 1,1,0.95))
    # lines(rep(x3, each=2), bracketpos*c(0.95, 1,1,0.95))
    par(xpd=F)  
    mtext("(B)", side=3, line=-0.25, adj=-0.15, cex=1, font=2)
    
    
    #===============================================================================
    # Suitability Persistence relative to temperature change (predicted baseline)
    #===============================================================================
    library(msir)
    hist.periods <- c("1991-2000", "1991-2017", "2001-2010", "2001-\n2017","2011-2017", "2017")
    ColScheme=c(2,1,3)
    
    plt=c(0.665, 0.995, 0.525, 0.99)
    ylim=c(0,1.1)
    par(mfrow=c(1,1), mar=c(3.25,3.25,0.1,0.1), mgp=c(1.5,0.25,0), plt = plt, new = TRUE)
    plot(0, xlim=c(0,7.9), ylim=ylim, yaxs="i", xaxs="i", col="white", xaxt="n", yaxt="n", 
         xlab=bquote(BC~mean~temperature~change~relative~to~"1961-90"~"("*degree*C*")"), 
         ylab="")
    axis(1, at=0:8, labels = 0:8, tck=0)
    axis(2, at=seq(0,1.2,0.2), labels = paste(seq(0,1.2,0.2)*100, "%", sep=""), las=2, tck=0)
    par(mgp=c(2,0.25,0), plt = plt, new = TRUE)
    title(ylab="BC mean suitability persistence             ")
    mtext("(C)", side=3, line=-0.75, adj= -0.1, cex=1, font=2)
    lines(c(-99,99), c(1,1), lwd=2, col="darkgrey")
    
    for(edatope.temp in edatopes){
      x <- c(0,MAT.change)
      y <- c(1,get(paste(metric, "mean.all", edatope.temp, sep=".")))
      l <- loess.sd(y~x, span=1, nsigma = 1)
      if(edatope.temp!=edatope) polygon(c(l$x, rev(l$x)), c(l$upper, rev(l$lower)), col=alpha(ColScheme[which(edatopes==edatope.temp)], 0.5), border=NA)
      l <- loess(y~x)
      par(xpd=T)
      text(max(x)-0.1, predict(l, max(x)), edatope.temp, pos=4, font= if(edatope.temp==edatope) 2 else 1, col=ColScheme[which(edatopes==edatope.temp)], cex= if(edatope.temp==edatope) 1.1 else 1)
      par(xpd=F)
    }
    
    temp.suit <- rep(NA, length(GCMs))
    temp.spp <- rep(NA, length(GCMs))
    for(rcp.temp in rcps){
      for(proj.year.temp in proj.years){
        temp.suit <- rep(NA, length(GCMs))
        temp.spp <- rep(NA, length(GCMs))
        
        x <- get(paste("MAT.change", rcp.temp, proj.year.temp, sep="."))
        y <- get(paste(metric, "mean", rcp.temp, proj.year.temp, edatope, sep=".")) 
        
        points(x,y, pch=c(21,22)[which(rcps==rcp.temp)], bg=c("black", "dodgerblue", "yellow")[which(proj.years==proj.year.temp)])
        # points(x,y2, col="red", pch=16)
        # points(0,mean(BGC.pred.ref!=BGC), pch=16, cex=1.3)
        # points(0,mean(zone.pred.ref!=zone), col="red", pch=16, cex=1.3)
      }
    }
    
    for(hist.year in hist.years[c(4)]){
      x <- get(paste("MAT.change", hist.year, sep="."))
      # y <- get(paste("SuitPersistence.mean", hist.year, edatope, sep="."))
      y <- get(paste(metric, "mean", hist.year, edatope, sep="."))
      points(x,y, cex=1.4, pch=2)
      text(x,y-0.01, paste(hist.periods[which(hist.years==hist.year)], ""), pos=1, cex=0.8, font=2)
      # print(y)
    }
    
    legend("bottomleft", legend=c("RCP4.5", "RCP8.5", "2011-2040", "2041-2070", "2071-2100"), pch=c(21,22, NA,NA,NA), 
           pt.bg=c("gray", "gray", NA,NA,NA), pt.cex=1.5, fill=c(NA, NA, "black", "dodgerblue", "yellow"), border = c(F, F, T, T, T), bty="n")
    
    # boxplot for focal period
    x <- c(0,MAT.change)
    x.focal <- MAT.change[which(seq.rcp==rcp & seq.proj.year==proj.year)]
    boxplot(x.focal, add=T, horizontal=TRUE, axes=FALSE, range=0, at=ylim[2]-0.05, boxwex = 0.12)
    text(max(x.focal), ylim[2]-0.05, paste(rcp.name[which(rcps==rcp)], ", ", proj.year.name[which(proj.years==proj.year)], sep=""), pos=4, cex=0.8)
    
    dev.off()
    print(edatope)
  }
  print(proj.year)
}





 
 
 
