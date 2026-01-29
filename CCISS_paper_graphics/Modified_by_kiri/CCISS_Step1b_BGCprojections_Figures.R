##======================================================================================
## Exploratory Data Analysis on the Biogeoclimatic Projections for the CCISS analysis
##======================================================================================

# Colin Mahony
# c_mahony@alumni.ubc.ca
# 778-288-4008


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

setwd("C:\\Users\\mahonyc.stu\\Documents\\Masters\\Research\\Publications\\2019_CCISS")

grid <- "BC2kmGrid"
grid.data <- read.csv(paste("InputData\\", grid, ".csv", sep = ""))

GCMs <-  c("ACCESS1-0","CanESM2","CCSM4","CESM1-CAM5","CNRM-CM5","CSIRO-Mk3-6-0", "GFDL-CM3","GISS-E2R", "HadGEM2-ES", "INM-CM4", "IPSL-CM5A-MR", "MIROC-ESM", "MIROC5", "MPI-ESM-LR","MRI-CGCM3")
rcps <- c("rcp45", "rcp85")
proj.years <- c(2025, 2055, 2085)
hist.years <- c(1995, 2004, 2005, 2009, 2014, 2017)
hist.year.name <- c("1991-2000", "1991-2017", "2001-2010", "2001-2017","2011-2017", "2017")

##Specify variables in model
model = "5.1"
VarList = c("AHM", "bFFP","CMD.total","DD5_sp","EMT","Eref_sm","EXT","FFP","MCMT","MSP",
            "PPT_JAS","PPT_MJ","PPT06","SHM","TD","Tmax_sp","Tmin_at","Tmin_sm","Tmin_wt",
            "PAS","CMD.def","CMDMax","eFFP","Eref09","MAT","PPT07","Tmin_sp")

BGCcolors <- read.csv("C:\\Users\\mahonyc.stu\\Documents\\Masters\\Research\\SpatialData\\BGCv10\\BGCv10\\BGCzone_Colorscheme.csv")

###Load random forest model
varset <- "16VAR"
fname=paste("Rcode\\BGCv11_AB_USA_LHC_", varset, "_SubZone_RFmodel.Rdata", sep="")
load(fname)


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

# # merge small units
# BGC[which(BGC=="CMAwh")] <- "CMAun"
# BGC[which(BGC=="MHun")] <- "MHunp"
# BGC[which(BGC=="ESSFxvw")] <- "ESSFxvp"
# BGC[which(BGC=="ESSFdcp")] <- "ESSFdcw"

#BGC zones
BGCcolors <- read.csv("C:\\Users\\mahonyc.stu\\Documents\\Masters\\Research\\SpatialData\\BGCv10\\BGCv10\\BGCzone_Colorscheme.csv")
zone <- rep(NA, length(BGC))
for(i in BGCcolors$zone){ zone[grep(i,BGC)] <- i }
table(zone)

#===============================================================================
# generic spatial data
#===============================================================================
### admin boundaries
bdy.bc <- readOGR("InputData\\BC_AB_US_Shp\\ProvincialOutline.shp")



#===============================================================================
# model predictions
#===============================================================================

## mapped BGC
points <- read.csv(paste("InputData\\",grid,".csv", sep=""))
BGC <- points$ID2
BGC <- gsub(" ","",BGC)  
zone <- rep(NA, length(BGC))
for(i in BGCcolors$zone){ zone[grep(i,BGC)] <- i }

## reference period BGC
BGC.pred.ref <- as.character(read.csv(paste("OutputData\\BGC.pred", grid, "ref.csv", sep="."))[,1])
zone.pred.ref <- rep(NA, length(BGC))
for(i in BGCcolors$zone){ zone.pred.ref[grep(i,BGC.pred.ref)] <- i }

# Historical BGC
for(hist.year in hist.years){
  BGC.pred <- as.character(read.csv(paste("OutputData\\BGC.pred", grid,hist.year,"csv", sep="."))[,1])
  assign(paste("BGC.pred", hist.year, sep="."), BGC.pred) #bgc projection
  print(hist.year)
}

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

# determine vote winner BGC and ensemble agreement (WARNING: takes about 2 minutes per rcp/proj.year)
# rcp=rcps[1]
for(rcp in rcps){
# proj.year=proj.years[2]
for(proj.year in proj.years){
temp <- as.data.frame(matrix(rep(NA, length(BGC.pred)*length(GCMs)), nrow=length(BGC.pred), ncol=length(GCMs)))
for(GCM in GCMs){
  BGC.pred <- get(paste("BGC.pred", GCM, rcp, proj.year, sep="."))
  #add votes to votes matrix
  temp[,which(GCMs==GCM)] <- BGC.pred
  # print(GCM)
}
vote.winner <- function(x){return(names(which(table(x)==max(table(x))))[1])}
agreement <- function(x){return(max(table(x)))}
assign(paste("BGC.pred.ensemble", rcp, proj.year, sep="."), apply(temp, 1, vote.winner))
# assign(paste("BGC.pred.agreement", rcp, proj.year, sep="."), apply(temp, 1, agreement))

    print(proj.year)
  }
  print(rcp)
}

# calculate percentage of ensemble projecting an exotic unit
# rcp=rcps[1]
for(rcp in rcps){
# proj.year=proj.years[2]
for(proj.year in proj.years){
exotic <- rep(0, length(BGC.pred.ref))
for(GCM in GCMs){
  BGC.pred <- get(paste("BGC.pred", GCM, rcp, proj.year, sep="."))
  temp <- rep(0, length(BGC.pred.ref))
  temp[-which(BGC.pred%in%unique(BGC))] <- 1
  exotic <- apply(cbind(exotic, temp), 1, sum)
  # print(GCM)
}
assign(paste("BGC.pred.exotic", rcp, proj.year, sep="."), exotic/length(GCMs))
    print(proj.year)
  }
  print(rcp)
}


#===============================================================================
#
# PLOTS
#
#===============================================================================

X <- dem

#===============================================================================
# BGC projections (Majority vote)
#===============================================================================

zones <- c("BG", "BWBS", "CDF", "CWH", "ESSF", "ICH", "IDF", "MH", "MS", "PP", "SBPS", "SBS", "SWB" )
ColScheme <- as.character(BGCcolors[match(zones,BGCcolors[,1]),5])

png(filename=paste("Results\\CCISS.BGCEDA.BGCprojections.vert","png",sep="."), type="cairo", units="in", width=6.5, height=8.5, pointsize=11, res=600)
par(mar=c(0.1,0.1, 0.1,0.1), mgp=c(2,0.25,0), mfrow=c(3,2))
  
# png(filename=paste("Results\\CCISS.BGCEDA.BGCprojections.horiz","png",sep="."), type="cairo", units="in", width=11, height=6.5, pointsize=11, res=600)
# par(mar=c(0.1,0.1, 0.1,0.1), mgp=c(2,0.25,0), mfrow=c(2,3))

#Mapped BGC zones
zone <- rep(NA, length(BGC))
for(i in BGCcolors$zone){ zone[grep(i,BGC)] <- i }
table(zone)

pred <- zone
pred[grep("IMA|CMA|BAFA", pred)] <- NA
values(X) <- NA
values(X) <- as.numeric(factor(pred))[plotOrder]
plot(X, xaxt="n", yaxt="n", col=ColScheme, legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE) 
plot(bdy.bc, add=T, lwd=0.4)
mtext("(a) BGCv10 Map", side=1, line=-1.5, adj=0.05, cex=1, font=2)
legend("left", legend=c(zones, "Alpine", "Exotic"), fill=c(ColScheme, "white", "black"), bty="n", cex=1.1, inset=-0.01)
# box()

#predicted BGC zones (ref period)
zone.pred.ref <- rep(NA, length(BGC.pred.ref))
for(i in BGCcolors$zone){ zone.pred.ref[grep(i,BGC.pred.ref)] <- i }
zone.pred.ref[-which(BGC.pred.ref%in%unique(BGC))] <- "Exotic"
pred <- zone.pred.ref
pred[grep("IMA|CMA|BAFA", pred)] <- NA
values(X) <- NA
values(X) <- as.numeric(factor(pred, levels=c(zones, "Exotic")))[plotOrder]
plot(X, xaxt="n", yaxt="n", col=c(ColScheme, "black"), legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE)
plot(bdy.bc, add=T, lwd=0.4)
mtext("(b) Predicted (1970s)", side=1, line=-1.5, adj=0.05, cex=1, font=2)

#predicted BGC zones (projected period)
for(j in 1:4){
  rcp=rcps[c(1,1,1,2)][j]
  proj.year=proj.years[c(1,2,3,3)][j]

  BGC.pred <- get(paste("BGC.pred.ensemble", rcp, proj.year, sep="."))
zone.pred <- rep(NA, length(BGC.pred))
for(i in BGCcolors$zone){ zone.pred[grep(i,BGC.pred)] <- i }
zone.pred[-which(BGC.pred%in%unique(BGC))] <- "Exotic"
pred <- zone.pred
pred[grep("IMA|CMA|BAFA", pred)] <- NA
values(X) <- NA
values(X) <- as.numeric(factor(pred, levels=c(zones, "Exotic")))[plotOrder]
plot(X, xaxt="n", yaxt="n", col=c(ColScheme, "black"), legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE)
plot(bdy.bc, add=T, lwd=0.4)
mtext(paste("(", letters[j+2], ") ", rcp, " ",proj.year, sep=""), side=1, line=-1.5, adj=0.05, cex=1, font=2)

}
dev.off()



#===============================================================================
# BGC projections (loop of all models)
#===============================================================================

############## get the full list of BGC zones
setwd("C:\\Users\\mahonyc.stu\\Documents\\Masters\\Research\\Publications\\2019_CCISS")

## parameters
grid <- "WNA2"
BGC.pred.ref.WNA <- read.csv(paste("OutputData\\BGC.pred.ref", grid, "csv", sep="."))[,1]
rcp="rcp45"
proj.year=2055
BGC.pred.WNA <- read.csv(paste("OutputData\\BGC.pred", grid, rcp, proj.year, "csv", sep="."))[,1]


zone.pred.ref.WNA <- gsub("[:a-z:]","",BGC.pred.ref.WNA) 
zone.pred.ref.WNA <- gsub("[:1-9:]","",zone.pred.ref.WNA) 

zone.pred.WNA <- gsub("[:a-z:]","",BGC.pred.WNA) 
zone.pred.WNA <- gsub("[:1-9:]","",zone.pred.WNA) 

zone <- sort(unique(c(zone.pred.ref.WNA, zone.pred.WNA)))
zone <- factor(zone, levels=zone)

#BGC zones
BGCcolors <- read.csv("C:\\Users\\mahonyc.stu\\Documents\\Masters\\Research\\SpatialData\\BGCv10\\BGCv10\\BGCzone_Colorscheme.csv")
colors = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)][-1]

ColScheme.zone <- rep(NA, length(zone))
ColScheme.zone[which(zone%in%BGCcolors$zone)] <- as.character(BGCcolors$HEX[match(zone[which(zone%in%BGCcolors$zone)], BGCcolors$zone)])
set.seed(8)
ColScheme.zone[-which(zone%in%BGCcolors$zone)] <- sample(colors, length(zone[-which(zone%in%BGCcolors$zone)]))
ColScheme.zone <- factor(ColScheme.zone, levels=ColScheme.zone)

################# 4-panel maps of all time periods for one model. 
GCM=GCMs[1]
for(GCM in GCMs){

png(filename=paste("Results\\CCISS.BGCEDA.BGCprojections", GCM,"png",sep="."), type="cairo", units="in", width=8.5, height=7.5, pointsize=11, res=600)
par(mar=c(0.1,0.1, 0.1,0.1), mgp=c(2,0.25,0), mfrow=c(2,2))

#predicted BGC zones (projected period)
for(j in 1:4){
  rcp=rcps[c(1,1,1,2)][j]
  proj.year=proj.years[c(1,2,3,3)][j]
  
  BGC.pred <- get(paste("BGC.pred", GCM, rcp, proj.year, sep="."))
  zone.pred <- rep(NA, length(BGC.pred))
  for(i in zone){ zone.pred[grep(i,BGC.pred)] <- i }
  exotic <- table(zone.pred[-which(zone.pred%in%BGCcolors$zone)])
  exotic <- exotic[exotic>20]
  exotic <- exotic[rev(order(exotic))]
  as.numeric(formatC(signif(exotic/length(zone.pred)*100,digits=3), digits=3,format="fg", flag="#"))
  pred <- zone.pred
  values(X) <- NA
  values(X) <- factor(pred, levels=zone)[plotOrder]
  values(X)[1:length(zone)] <- 1:length(zone) # this is a patch that is necessary to get the color scheme right.
  
  plot(X, xaxt="n", yaxt="n", col=alpha(ColScheme.zone, 1), legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE)
  plot(bdy.bc, add=T, lwd=0.4)
  mtext(paste("(", letters[j], ") ", rcp, " ",proj.year, sep=""), side=1, line=-4.5, adj=0.12, cex=1, font=2)
  mtext(GCM, side=1, line=-5.5, adj=0.115, cex=1.2, font=2)
  
  exotic.pct <- round(as.numeric(formatC(signif(exotic/length(zone.pred)*100,digits=3), digits=3,format="fg", flag="#")),2)
  legend("topright", legend=paste(names(exotic), " (", exotic.pct, "%)", sep=""), fill=alpha(ColScheme.zone[as.numeric(factor(names(exotic), zone))], 1), bty="n")
  

}
dev.off()
print(GCM)
}

################# 15-panel maps of all models for one time period
GCM=GCMs[1]
rcp=rcps[1]
proj.year=proj.years[2]

png(filename=paste("Results\\CCISS.BGCEDA.BGCprojections", rcp, proj.year,"png",sep="."), type="cairo", units="in", width=12.75, height=6.5, pointsize=10, res=600)
par(mar=c(0.1,0.1, 0.1,0.1), mgp=c(2,0.25,0), mfrow=c(3,5))

for(GCM in GCMs){
  
    BGC.pred <- get(paste("BGC.pred", GCM, rcp, proj.year, sep="."))
    zone.pred <- rep(NA, length(BGC.pred))
    for(i in zone){ zone.pred[grep(i,BGC.pred)] <- i }
    exotic <- table(zone.pred[-which(zone.pred%in%BGCcolors$zone)])
    exotic <- exotic[exotic>20]
    exotic <- exotic[rev(order(exotic))]
    as.numeric(formatC(signif(exotic/length(zone.pred)*100,digits=3), digits=3,format="fg", flag="#"))
    pred <- zone.pred
    values(X) <- NA
    values(X) <- factor(pred, levels=zone)[plotOrder]
    values(X)[1:length(zone)] <- 1:length(zone) # this is a patch that is necessary to get the color scheme right.
    
    plot(X, xaxt="n", yaxt="n", col=alpha(ColScheme.zone, 1), legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE)
    plot(bdy.bc, add=T, lwd=0.4)
    mtext(GCM, side=1, line=-1.5, adj=0.1, cex=1, font=2)
    
    exotic.pct <- round(as.numeric(formatC(signif(exotic/length(zone.pred)*100,digits=3), digits=3,format="fg", flag="#")),2)
    legend("topright", legend=paste(names(exotic), " (", exotic.pct, "%)", sep=""), fill=alpha(ColScheme.zone[as.numeric(factor(names(exotic), zone))], 1), bty="n")
    
    
  print(GCM)
}

dev.off()


################# 4-panel maps of historical time periods

  png(filename=paste("Results\\CCISS.BGCEDA.BGCprojections.histYears","png",sep="."), type="cairo", units="in", width=8.5, height=7.5, pointsize=11, res=600)
  par(mar=c(0.1,0.1, 0.1,0.1), mgp=c(2,0.25,0), mfrow=c(2,2))
  
  #predicted BGC zones 
  for(hist.year in hist.years[c(2,4,5,6)]){
   
    BGC.pred <- get(paste("BGC.pred", hist.year, sep="."))
    zone.pred <- rep(NA, length(BGC.pred))
    for(i in zone){ zone.pred[grep(i,BGC.pred)] <- i }
    exotic <- table(zone.pred[-which(zone.pred%in%BGCcolors$zone)])
    exotic <- exotic[exotic>20]
    exotic <- exotic[rev(order(exotic))]
    as.numeric(formatC(signif(exotic/length(zone.pred)*100,digits=3), digits=3,format="fg", flag="#"))
    pred <- zone.pred
    values(X) <- NA
    values(X) <- factor(pred, levels=zone)[plotOrder]
    values(X)[1:length(zone)] <- 1:length(zone) # this is a patch that is necessary to get the color scheme right.
    
    plot(X, xaxt="n", yaxt="n", col=alpha(ColScheme.zone, 1), legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE)
    plot(bdy.bc, add=T, lwd=0.4)
    mtext(paste("(", letters[j], ") ", rcp, " ",hist.year.name[which(hist.years==hist.year)], sep=""), side=1, line=-1.5, adj=0.12, cex=1, font=2)

    exotic.pct <- round(as.numeric(formatC(signif(exotic/length(zone.pred)*100,digits=3), digits=3,format="fg", flag="#")),2)
    if(length(exotic)>0) legend("topright", legend=paste(names(exotic), " (", exotic.pct, "%)", sep=""), fill=alpha(ColScheme.zone[as.numeric(factor(names(exotic), zone))], 1), bty="n")
    
    print(hist.year)
  }
  dev.off()


  
  #===============================================================================
  # Detail map of BGC unit projections (loop of all models)
  #===============================================================================
  
  ############## get the full list of BGC zones
  setwd("C:\\Users\\mahonyc.stu\\Documents\\Masters\\Research\\Publications\\2019_CCISS")
  
  ## parameters
  grid <- "WNA2"
  BGC.pred.ref.WNA <- read.csv(paste("OutputData\\BGC.pred.ref", grid, "csv", sep="."))[,1]
  rcp="rcp45"
  proj.year=2055
  BGC.pred.WNA <- read.csv(paste("OutputData\\BGC.pred", grid, rcp, proj.year, "csv", sep="."))[,1]
  
  BGC.WNA <- sort(unique(c(as.character(BGC.pred.ref.WNA), as.character(BGC.pred.WNA))))

  #BGC zones
  colors = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)][-1]
  
  ColScheme <- sample(colors, length(BGC.WNA))
  ColScheme <- factor(ColScheme, levels=ColScheme)
  
  ################# 4-panel maps of all time periods for one model. 
  GCM="ensemble"
  for(rcp in rcps){
    for(proj.year in proj.years){
  # for(hist.year in hist.years){
  # for(GCM in c("ensemble", GCMs)){
    
    png(filename=paste("Results\\CCISS.BGCEDA.BGCprojections.BGCUnits", GCM, rcp, proj.year,"png",sep="."), type="cairo", units="in", width=8.5, height=7.5, pointsize=8, res=600)
    # png(filename=paste("Results\\CCISS.BGCEDA.BGCprojections.BGCUnits.Reference.png",sep="."), type="cairo", units="in", width=8.5, height=7.5, pointsize=8, res=600)
    # png(filename=paste("Results\\CCISS.BGCEDA.BGCprojections.BGCUnits.Reference", hist.year,"png",sep="."), type="cairo", units="in", width=8.5, height=7.5, pointsize=8, res=600)
    par(mar=c(0.1,0.1, 0.1,0.1), mgp=c(2,0.25,0), mfrow=c(1,1))
    
    #predicted BGC zones (projected period)
    
    pred <- get(paste("BGC.pred", GCM, rcp, proj.year, sep="."))
    # pred <- BGC.pred.ref
    # pred <- get(paste("BGC.pred", hist.year, sep="."))
    values(X) <- NA
    values(X) <- factor(pred, levels=BGC.WNA)[plotOrder]
    values(X)[1:length(BGC.WNA)] <- 1:length(BGC.WNA) # this is a patch that is necessary to get the color scheme right.
    
    plot(X, xaxt="n", yaxt="n", col=alpha(ColScheme, 1), legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE)
    plot(bdy.bc, add=T, lwd=0.4)
    mtext(paste(GCM, rcp, proj.year), side=1, line=-5.5, adj=0.115, cex=1.2, font=2)
    # mtext(hist.year.name[which(hist.years==hist.year)], side=1, line=-5.5, adj=0.115, cex=1.2, font=2)
    # mtext("Reference period prediction", side=1, line=-5.5, adj=0.115, cex=1.2, font=2)
    
    bgcs <- BGC.WNA
    for(bgc in bgcs){
      pts <- which(BGC.WNA[values(X)]==bgc)
      if(length(pts)>1){ qs <- if(length(pts)>1000) c(0.05, 0.25, 0.5, 0.75, 0.95) else if(length(pts)>200) c(0.25, 0.5, 0.75) else c(0.5)
      for(q in qs){
        pt <- xyFromCell(X, pts[min(which(pts >= quantile(pts, q)))])
        points(pt, pch=21, bg=alpha(ColScheme[which(bgcs==bgc)], 1), cex=0.5, lwd=0.5)
        text(pt-c(10000, 0), bgc, pos=4, cex=0.25, font=2)
        # print(q)
      }
      }
      # print(paste(which(bgcs==bgc), "-", bgc))
    }
    box() 
    
    dev.off()
    # print(GCM)
  }
  }
  
  
  
#===============================================================================
# which models best match recent climate change? 
#===============================================================================

GCM <- "ensemble"
rcp <- rcps[1]
proj.year <- proj.years[1]
hist.year <- 2009

ctCompare.BGC <- vector()
ctCompare.zone <- vector()
for(GCM in GCMs){
  BGC.hist <- get(paste("BGC.pred", hist.year, sep="."))
  BGC.proj <- get(paste("BGC.pred", GCM, rcp, proj.year, sep="."))
  zone.hist <- rep(NA, length(BGC.hist))
  zone.proj <- rep(NA, length(BGC.proj))
  for(i in zone){ zone.hist[grep(i,BGC.hist)] <- i }
  for(i in zone){ zone.proj[grep(i,BGC.proj)] <- i }
  
  ct.BGC <- table(group=factor(BGC.hist, levels=BGC.WNA),class=factor(BGC.proj, levels=BGC.WNA))
  ct.zone <- table(group=factor(zone.hist, levels=zone),class=factor(zone.proj, levels=zone))
  
  ctCompare.BGC[which(GCMs==GCM)] <- sum(diag(ct.BGC))/sum(ct.BGC)
  ctCompare.zone[which(GCMs==GCM)] <- sum(diag(ct.zone))/sum(ct.zone)
print(GCM)
  }
  
png(filename=paste("Results\\CCISS.BGCEDA.RecentVsProjected","png",sep="."), type="cairo", units="in", width=6.5, height=6, pointsize=11, res=600)
par(mar=c(7,3,1,1), mgp=c(1.75, 0.25, 0))
barplot(t(cbind(ctCompare.BGC, ctCompare.zone)), col=c("grey80", "grey30"),  names.arg = GCMs, beside=T, las=2, tck=0, ylab="% of BC with same unit predicted", ylim=c(0,1))
legend("topleft",bty="n", fill=c("grey80", "grey30"), inset=0.05, legend=c("BGC subzone-variants", "BGC zones"), title = "Comparison of GCM projections (RCP4.5, 2011-2040)\nto recent BGC shifts (2001-2017)")
box()
dev.off()
#===============================================================================
# Ensemble Agreement
#===============================================================================

png(filename=paste("Results\\CCISS.BGCEDA.EnsembleAgreement","png",sep="."), type="cairo", units="in", width=6.5, height=6, pointsize=11, res=600)
par(mar=c(0.1,0.1, 0.1,0.1), mgp=c(2,0.25,0), mfrow=c(2,2))

for(j in 1:4){
  rcp=rcps[c(1,1,1,2)][j]
  proj.year=proj.years[c(1,2,3,3)][j]
  
  exotic <- get(paste("BGC.pred.agreement", rcp, proj.year, sep="."))
  
  values(X) <- NA
  values(X) <- exotic[plotOrder]
  breakpoints <- seq(0,1,1/length(GCMs)); length(breakpoints)
  # ColScheme <- c("white", colorRampPalette(brewer.pal(9,"YlOrBr"))(length(breakpoints)-2)); length(ColScheme)
  ColScheme <- colorRampPalette(c("white","gray90", "gray50", "gray10", "black"))(length(breakpoints)-1)
  
  plot(X, xaxt="n", yaxt="n", col=ColScheme, legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE) 
  box(col="white", lwd=1.5)
  plot(bdy.bc, add=T, lwd=0.4)
  if(j==1){
    xl <- 1600000; yb <- 1000000; xr <- 1700000; yt <- 1700000
    rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme, border = NA)
    rect(xl,  yb,  xr,  yt)
    text(rep(xr,length(breakpoints))[seq(1, length(breakpoints), 3)],seq(yb,yt,(yt-yb)/(length(breakpoints)-1))[seq(1, length(breakpoints), 3)],paste(round(breakpoints*100),"%", sep="")[seq(1, length(breakpoints), 3)],pos=4,cex=0.8,font=1)
    text(xl-40000, mean(c(yb,yt))-30000, paste("% of GCMs predicting\nthe majority unit"), srt=90, pos=3, cex=1, font=2)
  }
  mtext(paste(rcp, proj.year, sep="\n"), side=1, line=-3.5, adj=0.05, cex=1, font=2)
}
dev.off()

#===============================================================================
# exotic units in the ensemble projection
#===============================================================================

png(filename=paste("Results\\CCISS.BGCEDA.Exotic","png",sep="."), type="cairo", units="in", width=6.5, height=6, pointsize=11, res=600)
par(mar=c(0.1,0.1, 0.1,0.1), mgp=c(2,0.25,0), mfrow=c(2,2))

for(j in 1:4){
  rcp=rcps[c(1,1,1,2)][j]
  proj.year=proj.years[c(1,2,3,3)][j]
  
  exotic <- get(paste("BGC.pred.exotic", rcp, proj.year, sep="."))

  values(X) <- NA
  values(X) <- exotic[plotOrder]
  breakpoints <- seq(0,1,1/length(GCMs)); length(breakpoints)
  ColScheme <- c("white", colorRampPalette(brewer.pal(9,"YlOrBr"))(length(breakpoints)-2)); length(ColScheme)
  
  plot(X, xaxt="n", yaxt="n", col=ColScheme, legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE) 
  box(col="white", lwd=1.5)
  plot(bdy.bc, add=T, lwd=0.4)
  if(j==1){
    xl <- 1600000; yb <- 1000000; xr <- 1700000; yt <- 1700000
    rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme, border = NA)
    rect(xl,  yb,  xr,  yt)
    text(rep(xr,length(breakpoints))[seq(1, length(breakpoints), 3)],seq(yb,yt,(yt-yb)/(length(breakpoints)-1))[seq(1, length(breakpoints), 3)],paste(round(breakpoints*100),"%", sep="")[seq(1, length(breakpoints), 3)],pos=4,cex=0.8,font=1)
    text(xl-40000, mean(c(yb,yt))-30000, paste("% exotic units\nin ensemble projection"), srt=90, pos=3, cex=1, font=2)
  }
  mtext(paste(rcp, proj.year, sep="\n"), side=1, line=-3.5, adj=0.05, cex=1, font=2)
}
dev.off()

#===============================================================================
# reference period prediction error relative to unit size. 
#===============================================================================

par(mar=c(0.1,0.1, 0.1,0.1), mgp=c(2,0.75,0), mfrow=c(1,2))
values(X) <- NA
BGC.pred <- BGC.pred.ref
# BGC.pred <- get(paste("BGC.pred", GCM, rcp, proj.year, sep="."))
BGC.change <- (BGC.pred!=BGC)
values(X) <- BGC.change[plotOrder]
plot(X, col=c("grey", "black"), xaxt="n", yaxt="n", legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE) 
plot(bdy.bc, add=T, lwd=0.4)
# box(col="white")

BGC.size <- aggregate(BGC.change,by=list(BGC), FUN=length)[,-1]
BGC.sdMAT <- aggregate(MAT.ref,by=list(BGC), FUN=sd)[,-1]
BGC.pctError <- aggregate(BGC.change,by=list(BGC), FUN=mean)[,-1]
BGC.pctChange <- 1-table(BGC.pred)[match(names(table(BGC)), names(table(BGC.pred)))]/table(BGC)

png(filename=paste("Results\\CCISS.BGCEDA.PredError.SizeBias","png",sep="."), type="cairo", units="in", width=6.5, height=6.5, pointsize=12, res=600)
mat <- matrix(c(3,9,5,6,4,9,7,8,9,9,9,9,9,9,1,2),4, byrow=T)   #define the plotting order
layout(mat, widths=c(0.1,.1,1,1,1), heights=c(1,1,.05,.1))   #set up the multipanel plot

par(mar=c(0,0,0,0))


plot(1, type="n", axes=F, xlab="", ylab="")
text(1,1, "Mapped BGC unit size (# grid cells)", font=2,cex=1.2)

plot(1, type="n", axes=F, xlab="", ylab="")
text(1,1, "Standard deviation of BGC unit MAT range", font=2,cex=1.2)

plot(1, type="n", axes=F, xlab="", ylab="")
text(1,1, "Relative BGC unit area (Predicted/Mapped)", srt=90, font=2,cex=1.2)

plot(1, type="n", axes=F, xlab="", ylab="")
text(1,1, "RF model prediction error", srt=90, font=2,cex=1.2)

par(mar=c(0.1, 0.1, 0.1,0.1), mgp=c(2,0.25,0))

x <- log10(BGC.size)
y <- log(as.vector(unlist(BGC.pctChange+1)))
plot(x, y, col="white", yaxt="n", xaxt="n", ylab="", xlab="")
text(x, y, aggregate(BGC.change,by=list(BGC), FUN=length)[,1], cex=0.5)
par(xpd=T)
axis(2, at=log(c(0.25, 1/3, 0.5, 0.75, 1, 1.5, 2, 3, 4)), labels=c(0.25, 0.33, 0.5, 0.75, 1, 1.5, 2, 3, 4), tck=0, las=2)
par(xpd=F)
lines(c(-99,99), c(0,0), lty=2)

x <- log10(BGC.sdMAT)
y <- log(as.vector(unlist(BGC.pctChange+1)))
plot(x, y, col="white", yaxt="n", xaxt="n", ylab="", xlab="")
text(x, y, aggregate(BGC.change,by=list(BGC), FUN=length)[,1], cex=0.5)
lines(c(-99,99), c(0,0), lty=2)

x <- log10(BGC.size)
y <- as.vector(unlist(BGC.pctError))
plot(x, y, col="white", yaxt="n", xaxt="n", ylab="", xlab="")
text(x, y, aggregate(BGC.change,by=list(BGC), FUN=length)[,1], cex=0.5)
par(xpd=T)
axis(1, at=seq(1,10), labels=round(10^(seq(1,10))), tck=0)
axis(2, at=seq(0,1,0.1), labels=seq(0,1,0.1), tck=0, las=2)
par(xpd=F)

x <- log10(BGC.sdMAT)
y <- as.vector(unlist(BGC.pctError))
plot(x, y, col="white", yaxt="n", xaxt="n", ylab="", xlab="")
text(x, y, aggregate(BGC.change,by=list(BGC), FUN=length)[,1], cex=0.5)
par(xpd=T)
axis(1, at=seq(-10,10,0.2), labels=round(10^(seq(-10,10,0.2)),1))
par(xpd=F)

dev.off()

# #===============================================================================
# # size distribution of exotic units
# #===============================================================================
# 
# ## create a dem from the climateBC input data
# points.US <- read.csv(paste("InputData\\US2kmGrid.csv", sep=""))
# BGC.US <- points.US$id2
# BGC.US.sum <- table(BGC.US)
# BGC.US.sum <- BGC.US.sum[-which(names(BGC.US.sum)%in%BGC)]
# 
# points.AB <- read.csv(paste("InputData\\AB2kmGrid.csv", sep=""))
# BGC.AB <- points.AB$id2
# BGC.AB.sum <- table(BGC.AB)
# 
# png(filename=paste("Results\\CCISS.BGCEDA.JurisSizeDist","png",sep="."), type="cairo", units="in", width=6.5, height=4, pointsize=12, res=600)
# par(mfrow=c(1,1), mar=c(3,3,0.1,0.1), mgp=c(2,0.25,0))
# hist(log10(BGC.US.sum), freq=F, xaxt="n", xlab="Mapped BGC unit size (# grid cells)", border=NA, ylim=c(0,0.9), main="", tck=0, las=2, yaxs="i")
# polygon(density(log10(table(BGC))), col=alpha("grey", 0.5))
# polygon(density(log10(BGC.US.sum)), col=alpha("blue", 0.5))
# polygon(density(log10(BGC.AB.sum)), col=alpha("red", 0.5))
# axis(1, at=seq(1,10), labels=round(10^(seq(1,10))), tck=0)
# legend("topleft", legend=c("BC units", "AB units", "US units (predicted)"), fill=alpha(c("grey", "red", "blue"), 0.5), bty="n")
# box()
# dev.off()

#===============================================================================
# confusion tables for specified zones
#===============================================================================

rcp <- rcps[1]
proj.year <- proj.years[2]

BGC.pred <- get(paste("BGC.pred.ensemble", rcp, proj.year, sep="."))

group <- BGC
class <- BGC.pred

ct <- table(group=group,class=class)
ct.CWHinICH <- ct[grep("ICH", row.names(ct)),grep("CWH", colnames(ct))]
ct.CWHinICH <- ct.CWHinICH[which(apply(ct.CWHinICH, 1, sum)>0),]
ct.CWHinICH <- ct.CWHinICH[,which(apply(ct.CWHinICH, 2, sum)>0)]
write.csv(ct.CWHinICH, paste("OutputData//ct.CWHinICH", rcp, proj.year, "csv", sep="."))

ct.CWHinESSF <- ct[grep("ESSF", row.names(ct)),grep("CWH", colnames(ct))]
ct.CWHinESSF <- ct.CWHinESSF[which(apply(ct.CWHinESSF, 1, sum)>0),]
ct.CWHinESSF <- ct.CWHinESSF[,which(apply(ct.CWHinESSF, 2, sum)>0)]
sum.CWHinESSF <- apply(ct.CWHinESSF, 2, sum)
write.csv(ct.CWHinESSF, paste("OutputData//ct.CWHinESSF", rcp, proj.year, "csv", sep="."))

#===============================================================================
# climate space analysis
#===============================================================================

VarList <- row.names(importance(BGCmodel))


Columns = c("AHM", "bFFP",
            "CMD07","DD5_sp","EMT","Eref_sm","EXT","FFP","MCMT","MSP",
            "PPT07","PPT08", "PPT05","PPT06","PPT09", "SHM","TD","Tmax_sp","Tmin_at",
            "Tmin_sm","Tmin_wt", "PPT_at","PPT_wt", "PAS","eFFP",
            "Eref09","MAT","Tmin_sp","CMD")

# reference climatic means for BC units
grid <- "BC2kmGrid"
fplot=paste("InputData\\", grid, "_Normal_1961_1990MSY.csv", sep="")
Y0 <- fread(fplot, select = Columns, stringsAsFactors = FALSE, data.table = FALSE) #fread is faster than read.csv
#####generate some additional variables
Y0$PPT_MJ <- Y0$PPT05 + Y0$PPT06 # MaY/June precip
Y0$PPT_JAS <- Y0$PPT07 + Y0$PPT08 + Y0$PPT09 # July/Aug/Sept precip
Y0$PPT.dormant <- Y0$PPT_at + Y0$PPT_wt # for calculating spring deficit
Y0$CMD.def <- 500 - (Y0$PPT.dormant)# start of growing season deficit original value was 400 but 500 seems better
Y0$CMD.def [Y0$CMD.def < 0] <- 0 #negative values set to zero = no deficit
Y0$CMDMax <- Y0$CMD07
Y0$CMD.total <- Y0$CMD.def + Y0$CMD
Clim <- Y0
Clim.BGCs <- aggregate(Clim, by=list(BGC), FUN=mean, na.rm=T)
BGCs <- Clim.BGCs$Group.1
Clim.BGCs <- Clim.BGCs[,-1]
pca.BGCs <- prcomp(Clim.BGCs[,which(names(Clim.BGCs)%in%VarList)], scale=T)
pcscores.BGCs <- predict(pca.BGCs, Clim.BGCs)
zones <- rep(NA, length(BGCs))
for(i in BGCcolors$zone){ zones[grep(i,BGCs)] <- i }
zones <- factor(zones, BGCcolors$zone)

# reference climatic means for exotic units
grid <- "WNA2"
fplot=paste("InputData\\", grid, "_Normal_1961_1990MSY.csv", sep="")
Y0 <- fread(fplot, select = Columns, stringsAsFactors = FALSE, data.table = FALSE) #fread is faster than read.csv
#####generate some additional variables
Y0$PPT_MJ <- Y0$PPT05 + Y0$PPT06 # MaY/June precip
Y0$PPT_JAS <- Y0$PPT07 + Y0$PPT08 + Y0$PPT09 # July/Aug/Sept precip
Y0$PPT.dormant <- Y0$PPT_at + Y0$PPT_wt # for calculating spring deficit
Y0$CMD.def <- 500 - (Y0$PPT.dormant)# start of growing season deficit original value was 400 but 500 seems better
Y0$CMD.def [Y0$CMD.def < 0] <- 0 #negative values set to zero = no deficit
Y0$CMDMax <- Y0$CMD07
Y0$CMD.total <- Y0$CMD.def + Y0$CMD
Clim.WNA <- Y0
Clim.BGCs.WNA <- aggregate(Clim.WNA, by=list(BGC.pred.ref.WNA), FUN=mean, na.rm=T)
BGCs.WNA <- Clim.BGCs.WNA$Group.1
Clim.BGCs.WNA <- Clim.BGCs.WNA[,-1]
pca.BGCs.WNA <- prcomp(Clim.BGCs.WNA[-which(BGCs.WNA=="SASbo"),which(names(Clim.BGCs.WNA)%in%VarList)], scale=T) # SASbo has what seems to be an erroneously cold/wet climate
pcscores.BGCs.WNA <- predict(pca.BGCs.WNA, Clim.BGCs.WNA)
zones.WNA <- rep(NA, length(BGCs.WNA))
for(i in zone){ zones.WNA[grep(i,BGCs.WNA)] <- i }
zones.WNA <- factor(zones.WNA, zone)


# projected climatic means for bc units
grid <- "BC2kmGrid"
# for(rcp in rcps){
  for(GCM in GCMs){
    Y1 <- fread(paste("InputData\\", grid, "_", GCM, "_", rcp, "_BioVars.csv", sep=""), select = c("GCM", Columns), stringsAsFactors = FALSE, data.table = FALSE)
    Y1$PPT_MJ <- Y1$PPT05 + Y1$PPT06 # MaY/June precip
    Y1$PPT_JAS <- Y1$PPT07 + Y1$PPT08 + Y1$PPT09 # July/Aug/Sept precip
    Y1$PPT.dormant <- Y1$PPT_at + Y1$PPT_wt # for calculating spring deficit
    Y1$CMD.def <- 500 - (Y1$PPT.dormant)# start of growing season deficit original value was 400 but 500 seems better
    Y1$CMD.def [Y1$CMD.def < 0] <- 0 #negative values set to zero = no deficit
    Y1$CMDMax <- Y1$CMD07
    Y1$CMD.total <- Y1$CMD.def + Y1$CMD
    Clim <- Y1
    ## assign single vectors to RCPs and proj.years
    Ystr <- strsplit(Y1[,1], "_")
    Y4 <- matrix(unlist(Ystr), ncol=3, byrow=TRUE)
    # for(proj.year in proj.years){
    temp <- aggregate(Y1[which(Y4[,2]==rcp & Y4[,3]==proj.year),-1], by=list(BGC), FUN=mean, na.rm=T)[,-1]
    assign(paste("Clim.BGCs", GCM, rcp, proj.year, sep="."), temp)
      assign(paste("pcscores.BGCs", GCM, rcp, proj.year, sep="."), predict(pca.BGCs.WNA, temp))
      # }
    print(GCM)
  }
  # print(rcp)
# }

# #################################
# # Plot for selected raw variables
#   par(mar=c(3.25,3.25,0.1,0.1), mgp=c(2.25,0.25,0))
#   var1 <- "Tmin_sm"
# # for(var2 in VarList){
#   var2 <- "MCMT"
# log <- F
# bgc1 <- "CWHds1"
# bgc2 <- "ICHvk1"
# p1 <- which(BGCs==bgc1)
# p2 <- which(BGCs==bgc2)
# x <- Clim.BGCs[,which(names(Clim.BGCs)==var1)]
# y <- Clim.BGCs[,which(names(Clim.BGCs)==var2)]
# if(log==T) y <- log(y)
# 
# png(filename=paste("Results\\CCISS.BGCEDA.ClimateSpace", bgc1, bgc2, var1, var2, rcp, proj.year, "png",sep="."), type="cairo", units="in", width=6.5, height=6.5, pointsize=12, res=600)
# par(mar=c(3.25,3.25,0.1,0.1), mgp=c(2.25,0.25,0))
# plot(x,y, pch=21, bg=as.character(BGCcolors[match(zones,BGCcolors[,1]),5]), cex=1.3, tck=0, xlab=var1, ylab=var2)
# points(x[p1],y[p1], cex=2, pch=21, bg=as.character(BGCcolors[match(zones,BGCcolors[,1]),5])[p1], lwd=1.2)
# points(x[p2],y[p2], cex=2, pch=21, bg=as.character(BGCcolors[match(zones,BGCcolors[,1]),5])[p2], lwd=1.2)
# text(x[p2],y[p2], bgc2, pos=4, font=2, cex=0.8)
# for(GCM in GCMs){
#   x2 <- get(paste("Clim.BGCs", GCM, rcp, proj.year, sep="."))[,which(names(Clim.BGCs)==var1)]
#   y2 <- get(paste("Clim.BGCs", GCM, rcp, proj.year, sep="."))[,which(names(Clim.BGCs)==var2)]
#   if(log==T) y2 <- log(y2)
#   lines(c(x[p1],x2[p1]), c(y[p1],y2[p1]))
#   points(x2[p1], y2[p1], pch=16, cex=0.7)
# }
# points(x[p1],y[p1], cex=2, pch=21, bg=as.character(BGCcolors[match(zones,BGCcolors[,1]),5])[p1], lwd=1.2)
# text(x[p1],y[p1], bgc1, pos=2, font=2, cex=0.8)
# # }
# dev.off()

#################################
# Plot for PCs 1 and 2

PCs <- T
select <- which(BGCs.WNA%in%c("CWHds1", "ICHvk1", "IDFdk3", "SBSmc2"))

bgc1="CDFmm"
bgc2="CRFdhz"
bgc2="CRFdhz"
# for(bgc1 in BGCs.WNA[select]){
  
par(mar=c(3.25,3.25,0.1,0.1), mgp=c(2.25,0.25,0))
var1 <- if(PCs==T) 2 else "Tmin_sm"
var2 <- if(PCs==T) 3 else "MCMT"
log <- F
p1 <- which(BGCs.WNA==bgc1)
p2 <- which(BGCs.WNA==bgc2)
x <- if(PCs==T) pcscores.BGCs.WNA[,var1] else Clim.BGCs.WNA[,which(names(Clim.BGCs.WNA)==var1)]
y <- if(PCs==T) pcscores.BGCs.WNA[,var2] else Clim.BGCs.WNA[,which(names(Clim.BGCs.WNA)==var2)]
if(log==T) y <- log(y)

Clim.bgc1 <- if(bgc1 %in% BGCs) Clim[which(BGC==bgc1),] else Clim.WNA[which(BGC.pred.ref.WNA==bgc1),]
pcscores.bgc1 <- predict(pca.BGCs.WNA, Clim.bgc1)
x.bgc1 <- if(PCs==T) pcscores.bgc1[,var1] else Clim.bgc1[,which(names(Clim.WNA)==var1)]
y.bgc1 <- if(PCs==T) pcscores.bgc1[,var2] else Clim.bgc1[,which(names(Clim.WNA)==var2)]

Clim.bgc2 <- Clim.WNA[which(BGC.pred.ref.WNA==bgc2),]
pcscores.bgc2 <- predict(pca.BGCs.WNA, Clim.bgc2)
x.bgc2 <- pcscores.bgc2[,var1]
y.bgc2 <- pcscores.bgc2[,var2]

png(filename=paste("Results\\CCISS.BGCEDA.ClimateSpace.SpatialVar", varset, bgc1, var1,"X", var2, "png",sep="."), type="cairo", units="in", width=24, height=12, pointsize=10, res=300)
par(mar=c(3.25,3.25,0.1,0.1), mgp=c(2.25,0.25,0))
eqscplot(x[-which(BGCs.WNA=="SASbo")],y[-which(BGCs.WNA=="SASbo")], col="white", cex=2, tck=0, 
     xlab=if(PCs==T) paste("Climate PC", var1, sep="") else var1, 
     ylab=if(PCs==T) paste("Climate PC", var2, sep="") else var2)

points(x.bgc1, y.bgc1, pch=16, col=alpha(as.character(ColScheme.zone[match(zones.WNA,zone)])[p1],0.8))
points(x.bgc2, y.bgc2, pch=16, col=alpha(as.character(ColScheme.zone[match(zones.WNA,zone)])[p2],0.8))

# #spatial variation in bgc1
# z.ref <- kde2d(x.bgc1, y.bgc1, n=200, h=1)
# z.con <- contourLines(z.ref, levels=c(0.01))
# for(k in 1:length(z.con)) polygon(z.con[[k]]$x,z.con[[k]]$y, col=alpha(as.character(ColScheme.zone[match(zones.WNA,zone)])[p1],0.8), border=NA)
# 
# #spatial variation in bgc2
# z.ref <- kde2d(x.bgc2, y.bgc2, n=200, h=1)
# z.con <- contourLines(z.ref, levels=c(0.01))
# for(k in 1:length(z.con)) polygon(z.con[[k]]$x,z.con[[k]]$y, col=alpha(as.character(ColScheme.zone[match(zones.WNA,zone)])[p2],0.8), border=NA)

points(x,y, pch=21, bg=as.character(ColScheme.zone[match(zones.WNA,zone)]), cex=2)

text(x[-which(BGCs.WNA==bgc1)],y[-which(BGCs.WNA==bgc1)], BGCs.WNA[-which(BGCs.WNA==bgc1)], pos=4, font=2, cex=0.6)
# points(x[p2],y[p2], cex=2, pch=21, bg=as.character(ColScheme.zone[match(zones.WNA,zone)])[p2], lwd=1.2)
# text(x[p2],y[p2], bgc2, pos=4, font=1, cex=0.8)
for(GCM in GCMs){
  x2 <- if(PCs==T) get(paste("pcscores.BGCs", GCM, rcp, proj.year, sep="."))[,var1] else get(paste("Clim.BGCs", GCM, rcp, proj.year, sep="."))[,which(names(Clim.BGCs.WNA)==var1)]
  y2 <- if(PCs==T) get(paste("pcscores.BGCs", GCM, rcp, proj.year, sep="."))[,var2] else get(paste("Clim.BGCs", GCM, rcp, proj.year, sep="."))[,which(names(Clim.BGCs.WNA)==var2)]
  if(log==T) y2 <- log(y2)
  p1.bc <- which(BGCs==bgc1)
  lines(c(x[p1],x2[p1.bc]), c(y[p1],y2[p1.bc]))
  bgc.pred.p1 <- as.character(predict(BGCmodel, get(paste("Clim.BGCs", GCM, rcp, proj.year, sep="."))[p1.bc,]))
  zone.pred.p1 <- gsub("[:a-z:]","",bgc.pred.p1) 
  zone.pred.p1 <- gsub("[:1-9:]","",zone.pred.p1) 
  points(x2[p1.bc], y2[p1.bc], pch=22, cex=1.2, bg=alpha(as.character(ColScheme.zone[match(zone.pred.p1,zone)]),1))
  text(x2[p1.bc], y2[p1.bc], bgc.pred.p1, cex=0.7, pos=4, col="grey")
}
points(x[p1],y[p1], cex=3, pch=21, bg=as.character(ColScheme.zone[match(zones.WNA,zone)])[p1], lwd=1.2)
text(x[p1],y[p1], bgc1, pos=4, font=2, cex=1.2)
legend("topleft", legend=bgc1, cex=1.5, bty="n")
dev.off()

print(paste(which(BGCs.WNA==bgc1), bgc1))
# }

