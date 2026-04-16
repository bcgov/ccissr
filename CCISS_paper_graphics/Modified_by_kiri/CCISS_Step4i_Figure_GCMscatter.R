

##======================================================================================
## CCISS Publication Scripts
## Step 4i - Figure - plot of temp and precip change by model 
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
# GCMs <-  GCMs[1:2] #XXX for testing purposes
proj.year.name <- c("2021-2040", "2041-2070", "2071-2100")

#===============================================================================
# generate the vector of mapped BGCs
#===============================================================================

points <- read.csv(paste("InputData\\",grid,".csv", sep=""))
BGC <- points$ID2
BGC <- gsub(" ","",BGC)  

#BGC zones
BGCcolors <- read.csv("C:\\Colin\\Projects\\2019_CCISS\\InputData\\BGCzone_Colorscheme.csv")
zone <- rep(NA, length(BGC))
for(i in BGCcolors$zone){ zone[grep(i,BGC)] <- i }
table(zone)

region <- c("N", "C", "S", "N", "N", "S", NA, "C", NA, "S", "S", "S", "S", "C", NA, "C")
BGC.region <- region[match(zone, BGCcolors$zone)]
BGC.regions <- c("C", "N", "S")
table(zone, BGC.region)
#===============================================================================
# Climate Data for reference period
#===============================================================================


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
MSP.ref <- Y0$PPT_MJ+Y0$PPT_JAS
PPT_wt.ref <- Y0$PPT_wt
MAT.mean.ref <- mean(MAT.ref, na.rm=T)
MSP.mean.ref <- mean(MSP.ref, na.rm=T)
PPT_wt.mean.ref <- mean(PPT_wt.ref, na.rm=T)



#===============================================================================
# Climate Data for future period
#===============================================================================
rcp="rcp45"
GCM=GCMs[1]
rcps <- c("rcp45", "rcp85")

for(GCM in GCMs){
  Y1 <- fread(paste("InputData\\", grid, "_", GCM, "_BioVars.csv", sep=""), stringsAsFactors = FALSE, data.table = FALSE)
  
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
      assign(paste("MSP", GCM, rcp, proj.year, sep="."), Y1$PPT_MJ[which(Y4[,2]==rcp & Y4[,3]==proj.year)]+Y1$PPT_JAS[which(Y4[,2]==rcp & Y4[,3]==proj.year)])
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
      MAT.change <- MAT-MAT.ref
      MSP.change <- MSP/MSP.ref
      PPT_wt.change <- PPT_wt/PPT_wt.ref
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

# same thing for specific regions
for(region in BGC.regions){
  select <- which(BGC.region==region)
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
      MAT.change <- MAT[select]-MAT.ref[select]
      MSP.change <- MSP[select]/MSP.ref[select]
      PPT_wt.change <- PPT_wt[select]/PPT_wt.ref[select]
      MAT.change.mean[which(GCMs==GCM)] <- mean(MAT.change, na.rm=T)
      MSP.change.mean[which(GCMs==GCM)] <- mean(MSP.change, na.rm=T)
      PPT_wt.change.mean[which(GCMs==GCM)] <- mean(PPT_wt.change, na.rm=T)
    }
    assign(paste("MAT.change", rcp, proj.year, region, sep="."), MAT.change.mean)
    assign(paste("MSP.change", rcp, proj.year, region, sep="."), MSP.change.mean)
    assign(paste("PPT_wt.change", rcp, proj.year, region, sep="."), PPT_wt.change.mean)
    print(proj.year)
  }
  print(rcp)
}
}

#===============================================================================
# change in climate variables for the ensemble
#===============================================================================

png(filename=paste("Results\\GCMscatter\\CCISS.BGCEDA.TempPrecipChange","png",sep="."), type="cairo", units="in", width=6.5, height=4, pointsize=9, res=600)

mat <- matrix(c(13,13,10,11,12,1,13,4,6,8,2,13,5,7,9,13,13,3,3,3),4, byrow=T)   #define the plotting order
layout(mat, widths=c(0.1,.1,1,1,1), heights=c(0.1,1,1,.2))   #set up the multipanel plot

par(mar=c(0,0,0,0))

plot(1, type="n", axes=F, xlab="", ylab="")
text(1,1, "Change in summer precipitation", srt=90, font=2,cex=1.2)

plot(1, type="n", axes=F, xlab="", ylab="")
text(1,1, "Change in winter precipitation", srt=90, font=2,cex=1.2)

plot(1, type="n", axes=F, xlab="", ylab="")
text(1,1, bquote(bold(Change~"in"~mean~annual~temperature~"("*degree*C*")")), font=2,cex=1.2)

par(mar=c(0.1,0.1,0.1,0.1), mgp=c(2,0.25,0))
for(proj.year in proj.years){
  x <- get(paste("MAT.change", rcps[1], proj.year, sep="."))
  y <- get(paste("MSP.change", rcps[1], proj.year, sep="."))
  x2 <- get(paste("MAT.change", rcps[2], proj.year, sep="."))
  y2 <- get(paste("MSP.change", rcps[2], proj.year, sep="."))
  plot(x,y,col="white", tck=0, xlim=c(0,8), ylim=c(0.85,1.22), yaxt="n", xaxt="n", xlab="", ylab="")
  lines(c(-9,9), c(1,1), col="gray", lty=2)
  text(x,y,GCMs, cex=0.6)
  text(x2,y2,GCMs, cex=0.6, col="dodgerblue")
  par(xpd=T)
  if(proj.year==proj.years[1]){
    axis(2, at=seq(0,2,0.1), labels=paste(seq(0,2,0.1)*100-100, "%", sep=""), las=2, tck=0)
    legend("topright", legend=c("RCP4.5", "RCP8.5"), cex=1.4, text.col = c("black", "dodgerblue"), bty="n")
  }
  mtext(paste("(", LETTERS[which(proj.years==proj.year)], ")", sep=""), side=3, line=-1.5, adj=0.05, cex=1, font=2)
  par(xpd=F)
  
  
  x <- get(paste("MAT.change", rcps[1], proj.year, sep="."))
  y <- get(paste("PPT_wt.change", rcps[1], proj.year, sep="."))
  x2 <- get(paste("MAT.change", rcps[2], proj.year, sep="."))
  y2 <- get(paste("PPT_wt.change", rcps[2], proj.year, sep="."))
  plot(x,y,col="white", tck=0, xlim=c(0,8), ylim=c(0.95,1.34), yaxt="n", xaxt="n", xlab="", ylab="")
  lines(c(-9,9), c(1,1), col="gray", lty=2)
  text(x,y,GCMs, cex=0.6)
  text(x2,y2,GCMs, cex=0.6, col="dodgerblue")
  par(xpd=T)
  if(proj.year==proj.years[1]) axis(2, at=seq(0,2,0.1), labels=paste(seq(0,2,0.1)*100-100, "%", sep=""), las=2, tck=0)
  axis(1, at=seq(0,8,2), labels=seq(0,8,2), tck=0)
  mtext(paste("(", LETTERS[4:6][which(proj.years==proj.year)], ")", sep=""), side=3, line=-1.5, adj=0.05, cex=1, font=2)
  par(xpd=F)
}

for(proj.year in proj.years){
  plot(1, type="n", axes=F, xlab="", ylab="")
  text(1,1, proj.year.name[which(proj.years==proj.year)], font=2,cex=1.2)
}
dev.off()

#===============================================================================
# change in climate variables for the ensemble
#===============================================================================
region.names <- c("Coast", "North", "South")
for(region in BGC.regions){
  
png(filename=paste("Results\\GCMscatter\\CCISS.BGCEDA.TempPrecipChange", region.names[which(BGC.regions==region)],"png",sep="."), type="cairo", units="in", width=6.5, height=4, pointsize=9, res=600)

mat <- matrix(c(13,13,10,11,12,1,13,4,6,8,2,13,5,7,9,13,13,3,3,3),4, byrow=T)   #define the plotting order
layout(mat, widths=c(0.1,.1,1,1,1), heights=c(0.1,1,1,.2))   #set up the multipanel plot

par(mar=c(0,0,0,0))

plot(1, type="n", axes=F, xlab="", ylab="")
text(1,1, "Change in summer precipitation", srt=90, font=2,cex=1.2)

plot(1, type="n", axes=F, xlab="", ylab="")
text(1,1, "Change in winter precipitation", srt=90, font=2,cex=1.2)

plot(1, type="n", axes=F, xlab="", ylab="")
text(1,1, bquote(bold(Change~"in"~mean~annual~temperature~"("*degree*C*")")), font=2,cex=1.2)

par(mar=c(0.1,0.1,0.1,0.1), mgp=c(2,0.25,0))
for(proj.year in proj.years){
  x <- get(paste("MAT.change", rcps[1], proj.year, region, sep="."))
  y <- get(paste("MSP.change", rcps[1], proj.year, region, sep="."))
  x2 <- get(paste("MAT.change", rcps[2], proj.year, region, sep="."))
  y2 <- get(paste("MSP.change", rcps[2], proj.year, region, sep="."))
  plot(x,y,col="white", tck=0, xlim=c(0,8), ylim=c(0.85,1.20), yaxt="n", xaxt="n", xlab="", ylab="")
  lines(c(-9,9), c(1,1), col="gray", lty=2)
  text(x,y,GCMs, cex=0.6)
  text(x2,y2,GCMs, cex=0.6, col="dodgerblue")
  par(xpd=T)
  if(proj.year==proj.years[1]){
    axis(2, at=seq(0,2,0.1), labels=paste(seq(0,2,0.1)*100-100, "%", sep=""), las=2, tck=0)
    legend("topright", legend=c("RCP4.5", "RCP8.5"), cex=1.4, text.col = c("black", "dodgerblue"), bty="n")
  }
  mtext(paste("(", LETTERS[which(proj.years==proj.year)], ")", sep=""), side=3, line=-1.5, adj=0.05, cex=1, font=2)
  par(xpd=F)
  
  x <- get(paste("MAT.change", rcps[1], proj.year, region, sep="."))
  y <- get(paste("PPT_wt.change", rcps[1], proj.year, region, sep="."))
  x2 <- get(paste("MAT.change", rcps[2], proj.year, region, sep="."))
  y2 <- get(paste("PPT_wt.change", rcps[2], proj.year, region, sep="."))
  plot(x,y,col="white", tck=0, xlim=c(0,8), ylim=c(0.95,1.32), yaxt="n", xaxt="n", xlab="", ylab="")
  lines(c(-9,9), c(1,1), col="gray", lty=2)
  text(x,y,GCMs, cex=0.6)
  text(x2,y2,GCMs, cex=0.6, col="dodgerblue")
  par(xpd=T)
  if(proj.year==proj.years[1]) axis(2, at=seq(0,2,0.1), labels=paste(seq(0,2,0.1)*100-100, "%", sep=""), las=2, tck=0)
  axis(1, at=seq(0,8,2), labels=seq(0,8,2), tck=0)
  mtext(paste("(", LETTERS[4:6][which(proj.years==proj.year)], ")", sep=""), side=3, line=-1.5, adj=0.05, cex=1, font=2)
  par(xpd=F)
}

for(proj.year in proj.years){
  plot(1, type="n", axes=F, xlab="", ylab="")
  text(1,1, proj.year.name[which(proj.years==proj.year)], font=2,cex=1.2)
}
dev.off()

}
