
##======================================================================================
## CCISS Publication Scripts
## Step 4e - Figures - Recent periods: Suitability for individual species and groups of species
##======================================================================================

# Colin Mahony
# c_mahony@alumni.ubc.ca
# 778-288-4008
# July 21, 2019

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

GCMs <-  c("ACCESS1-0","CanESM2","CCSM4","CESM1-CAM5","CNRM-CM5","CSIRO-Mk3-6-0", "GFDL-CM3","GISS-E2R", "HadGEM2-ES", "INM-CM4", "IPSL-CM5A-MR", "MIROC-ESM", "MIROC5", "MPI-ESM-LR","MRI-CGCM3")    
rcps <- c("rcp45", "rcp85")
proj.years <- c(2025, 2055, 2085)
hist.years <- c(1995, 2004, 2005, 2009, 2014, 2017)
hist.year.name <- c("1991-2000", "1991-2017", "2001-2010", "2001-2017","2011-2017", "2017")
edatopes<- c("B2", "C4", "D6")

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

#===============================================================================
# Manuscript plot of C4 for four species
#===============================================================================
hist.year <- 2017
edatope="C4"
spps <- c("Pl", "Fd", "Cw", "Sx")
# spps <- c("Ba", "Bl", "Bg")
# spps <- c("Yc", "Pa", "Hm")
# spps <- c("Lw", "Hw", "Py")
# spps <- c("Dr", "Ep", "At")

# x11(width=6.5, height=6, pointsize=10)

png(filename=paste("Results\\Manu_Suitability_Groups\\Recent\\CCISS.manu.Suitability.4panel",spps[1],spps[2],spps[3],spps[4], edatope, hist.year,"png",sep="."), type="cairo", units="in", width=6.5, height=6, pointsize=10, res=600)
par(mar=c(0,0,0,0), mfrow=c(1,1), bg="white", xpd=T)

#=============================
## Base plot
par(plt=c(0,1,0,1), bg=NA)
plot(0, col="white", xaxt="n", yaxt="n", xlab="", ylab="")
box(col="white")

## plot position
x1 <- c(0.05, 0.45,0.05,0.45)
x2 <- c(0.6,1,0.6,1)
y1 <- c(0.5, 0.5,0,0)
y2 <- c(1,1,0.5,0.5)

#=============================
## loop for each species
for(spp in spps){
  par(plt = as.vector(rbind(x1, x2, y1, y2)[,which(spps==spp)]), new = TRUE)
  
  RefSuit <- read.csv(paste("OutputData\\Suit.ref", grid, spp, edatope, "csv", sep="."))[,1]
  outRange.base <- RefSuit==5
  RefSuit[RefSuit==5] <- 4
  RefSuit[is.na(RefSuit)] <- 4
  
  ProjSuit <- read.csv(paste("OutputData\\Suit", grid, hist.year, spp, edatope, "csv", sep="."))[,1]
  ProjSuit[ProjSuit==5] <- 4
  ProjSuit[is.na(ProjSuit)] <- 4
  
  changeSuit <- RefSuit-ProjSuit
  outRange <- outRange.base
  outRange[which(changeSuit!=0)] <- FALSE
  changeSuit[outRange==T] <- NA
  
  breakpoints <- seq(-3.5,3.5,1); length(breakpoints)
  labels <- c("-3", "no\nchange", "+3")
  ColScheme <- c(brewer.pal(11,"RdBu")[c(1,2,3)], "grey80", brewer.pal(11,"RdBu")[c(9,10,11)]); length(ColScheme)
  
  values(X) <- changeSuit[plotOrder]
  plot(bdy.bc, border="black", lwd=0.4)
  image(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks=breakpoints, maxpixels= ncell(X))
  plot(bdy.bc, add=T, border="black", lwd=0.4)
  if(spp==spps[1]){
  xl <- 300000; yb <- 300000; xr <- 400000; yt <- 1200000
  rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme)
  text(rep(xr,length(labels)),seq(yb,yt,(yt-yb)/(length(labels)*5-1))[c(2,8,14)],labels,pos=4,cex=1,font=1)
  text(xl-30000, mean(c(yb,yt))-30000, paste("Change in suitability\n(", hist.year.name[which(hist.years==hist.year)], ")", sep=""), srt=90, pos=3, cex=1, font=2)
  }
  mtext(paste("(", LETTERS[which(spps==spp)],") ", spp, sep=""), side=3, line=-3.75, adj=0.05, cex=1, font=2)
  
print(spp)
}
dev.off()

#===============================================================================
# Manuscript plot of all edatopes for one species
#===============================================================================
rcp=rcps[1]
hist.year=hist.years[4]

# select the species to run the analysis on
spps <- unique(S1$Spp)
spps <- spps[-which(spps=="X")]
spps.candidate <- spps.lookup$TreeCode[-which(spps.lookup$Exclude=="x")]
spps <- spps[which(spps%in%spps.candidate)] 

for(spp in spps){
png(filename=paste("Results\\Manu_Suitability_Indiv\\Recent\\CCISS.suitability.edatopes",spp, hist.year,"png",sep="."), type="cairo", units="in", width=6.5, height=8.5, pointsize=12, res=600)
par(mar=c(0,0,0,0), mfrow=c(3,1), bg="white")

for(edatope in edatopes){
  
  RefSuit <- read.csv(paste("OutputData\\Suit.ref", grid, spp, edatope, "csv", sep="."))[,1]
  outRange.base <- RefSuit==5
  RefSuit[RefSuit==5] <- 4
  RefSuit[is.na(RefSuit)] <- 4
  
  ProjSuit <- read.csv(paste("OutputData\\Suit", grid, hist.year, spp, edatope, "csv", sep="."))[,1]
  ProjSuit[ProjSuit==5] <- 4
  ProjSuit[is.na(ProjSuit)] <- 4
  
  changeSuit <- RefSuit-ProjSuit
  outRange <- outRange.base
  outRange[which(changeSuit!=0)] <- FALSE
  changeSuit[outRange==T] <- NA
  
  par(plt=c(0,1,0,1), bg=NA)
  plot(0, col="white", xaxt="n", yaxt="n", xlab="", ylab="")
  Common <- as.character(spps.lookup$EnglishName[which(spps.lookup$TreeCode==spp)])
  Latin <- as.character(spps.lookup$ScientificName[which(spps.lookup$TreeCode==spp)])
  panel <- paste("(", LETTERS[which(edatopes==edatope)],")", sep="")
  mtext(paste(panel, " Site type: ", edatope, " (", edatope.name[which(edatopes==edatope)], ")", sep=""), side=3, line=-1.75, adj=0.01, cex=0.8, font=2)
  mtext(if(spp%in%spps.lookup$TreeCode) bquote(bold(.(spp))~"-"~.(Common)) else bquote(bold(.(spp))), side=3, line=-3, adj=0.01, cex=0.7, font=1)
  box()
  
  breakpoints <- seq(-3.5,3.5,1); length(breakpoints)
  labels <- c("-3", "no\nchange", "+3")
  ColScheme <- c(brewer.pal(11,"RdBu")[c(1,2,3)], "grey80", brewer.pal(11,"RdBu")[c(9,10,11)]); length(ColScheme)
  
  par(plt = c(0.1,1,0,1), new = TRUE)
  values(X) <- changeSuit[plotOrder]
  plot(bdy.bc, border="black", lwd=0.4)
  image(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks=breakpoints, maxpixels= ncell(X))
  plot(bdy.bc, add=T, border="black", lwd=0.4)
  # if(spp==spps[1]){
  xl <- 300000; yb <- 1000000; xr <- 400000; yt <- 1550000
  rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme)
  text(rep(xr,length(labels)),seq(yb,yt,(yt-yb)/(length(labels)*5-1))[c(2,8,14)],labels,pos=4,cex=1,font=1)
  text(xl-30000, mean(c(yb,yt))-30000, paste("Projected change\nin suitability (", hist.year.name[which(hist.years==hist.year)], ")", sep=""), srt=90, pos=3, cex=1, font=2)
  # }
  
  
  breakseq <- c(0.5,1.5,2.5,3.5,5)
  ColScheme <- c(brewer.pal(9,"Greys")[9], brewer.pal(9,"Greens")[7], brewer.pal(9,"Greens")[4], "white")
  length(ColScheme)
  
  par(plt = c(0, 0.375, 0, 0.75), new = TRUE)
  values(X) <- RefSuit[plotOrder]
  plot(bdy.bc, border="black", lwd=0.4)
  image(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks=breakseq, maxpixels= ncell(X))
  plot(bdy.bc, add=T, border="black", lwd=0.4)
  # if(spp==spps[1]){
  legend("bottomleft", legend=c("1 (primary)", "2 (secondary)", "3 (tertiary)"), 
         fill=ColScheme, bty="n", cex=0.9, title="Historical suitability", inset=0.015)
  # }
  # box()
  
  # map of binary appearance/disappearance
  binary <- rep(0, length(changeSuit))
  binary[outRange.base==T] <- NA
  binary[outRange.base][ProjSuit[outRange.base]<4] <- 1
  binary[outRange.base==F][ProjSuit[outRange.base==F]==4] <- -1
  values(X) <- binary[plotOrder]

  labels <- c("Disappearing", "No Change", "Appearing")
  ColScheme <- c(brewer.pal(11,"RdBu")[2], "grey85", brewer.pal(11,"RdBu")[10]); length(ColScheme)
  
  par(plt = c(0.625, 1, 0.25, 1), new = TRUE)
  plot(bdy.bc, border="black", lwd=0.4)
  image(X, add=T, xaxt="n", yaxt="n", col=ColScheme, maxpixels= ncell(X))
  
  legend("topright", legend=labels, fill=ColScheme, bty="n", cex=0.9)

  # print(edatope)
}
dev.off()
print(spp)
}


#===============================================================================
# Basic plot for a single species
#===============================================================================
# spps.name <- c("Yellow cedar", "Interior spruce", "Ponderosa pine", "Lodgepole pine", "Western larch", "Western hemlock", "Mountain hemlock", "Douglas-fir", "Red alder", "Western redcedar", "Subalpine fir", "Grand fir")
spps.lookup <- read.csv("InputData\\Tree speciesand codes_2.0.csv")
edatope.name <- c("Subxeric-poor", "Mesic-medium", "Hygric-rich")
rcp=rcps[1]
hist.year=hist.years[4]

edatope="C4"
# for(edatope in edatopes){

spp=spps[1]
for(spp in spps){
  
  png(filename=paste("Results\\CCISS.Suitability",spp, edatope, hist.year,"png",sep="."), type="cairo", units="in", width=8, height=2.35, pointsize=10, res=600)
  par(mar=c(0,0,0,0), mfrow=c(1,3), bg="white")
  
  RefSuit <- read.csv(paste("OutputData\\SI_Suitability_Basic\\Recent\\Suit.ref", grid, spp, edatope, "csv", sep="."))[,1]
  outRange.base <- RefSuit==5
  RefSuit[RefSuit==5] <- 4
  RefSuit[is.na(RefSuit)] <- 4
  
  # suitability for hist year
  ProjSuit <- read.csv(paste("OutputData\\Suit", grid, hist.year, spp, edatope, "csv", sep="."))[,1]
  ProjSuit[ProjSuit==5] <- 4
  ProjSuit[is.na(ProjSuit)] <- 4

  changeSuit <- RefSuit-ProjSuit
  outRange <- outRange.base
  outRange[which(changeSuit!=0)] <- FALSE
  changeSuit[outRange==T] <- NA
  
  breakseq <- c(0.5,1.5,2.5,3.5,5)
  ColScheme <- c(brewer.pal(9,"Greys")[9], brewer.pal(9,"Greens")[7], brewer.pal(9,"Greens")[4], "white")
  length(ColScheme)
  
  values(X) <- RefSuit[plotOrder]
  plot(bdy.bc, border="black", lwd=0.4)
  image(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks=breakseq, maxpixels= ncell(X))
  plot(bdy.bc, add=T, border="black", lwd=0.4)
  # if(spp==spps[1]){
  legend("topright", legend=c("1 (primary)", "2 (secondary)", "3 (tertiary)"), 
         fill=ColScheme, bty="n", cex=0.9, title="Historical suitability", inset=0.015)
  # }
  # box()
  Common <- as.character(spps.lookup$EnglishName[which(spps.lookup$TreeCode==spp)])
  Latin <- as.character(spps.lookup$ScientificName[which(spps.lookup$TreeCode==spp)])
  panel <- paste("(", LETTERS[which(spps==spp)],")", sep="")
  mtext(if(spp%in%spps.lookup$TreeCode) bquote(bold(.(spp))~"-"~.(Common)) else bquote(.(panel)~bold(.(spp))),
        side=1, line=-3, adj=0.01, cex=0.8, font=2)
  # mtext(if(spp%in%spps.lookup$TreeCode) bquote(.(panel)~bold(.(spp))~"-"~.(Common)~"("*italic(.(Latin)*")")) else bquote(.(panel)~bold(.(spp))),
  #       side=3, line=-1.75, adj=0.01, cex=0.8, font=2)
  mtext(paste("Site type: ", edatope, " (", edatope.name[which(edatopes==edatope)], ")", sep=""), side=1, line=-1.5, adj=0.01, cex=0.7, font=1)
  
  
  breakpoints <- seq(-3.5,3.5,1); length(breakpoints)
  labels <- c("-3", "no\nchange", "+3")
  ColScheme <- c(brewer.pal(11,"RdBu")[c(1,2,3)], "grey80", brewer.pal(11,"RdBu")[c(9,10,11)]); length(ColScheme)
  
  values(X) <- changeSuit[plotOrder]
  plot(bdy.bc, border="black", lwd=0.4)
  image(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks=breakpoints, maxpixels= ncell(X))
  plot(bdy.bc, add=T, border="black", lwd=0.4)
  # if(spp==spps[1]){
  xl <- 1600000; yb <- 1000000; xr <- 1700000; yt <- 1700000
  rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme)
  text(rep(xr,length(labels)),seq(yb,yt,(yt-yb)/(length(labels)*5-1))[c(2,8,14)],labels,pos=4,cex=1,font=1)
  text(xl-30000, mean(c(yb,yt))-30000, paste("Projected change\nin suitability (", hist.year.name[which(hist.years==hist.year)], ")", sep=""), srt=90, pos=3, cex=1, font=2)
  # }
  
  # map of binary appearance/disappearance
  binary <- rep(0, length(changeSuit))
  binary[outRange.base==T] <- NA
  binary[outRange.base][ProjSuit[outRange.base]<4] <- 1
  binary[outRange.base==F][ProjSuit[outRange.base==F]==4] <- -1
  values(X) <- binary[plotOrder]
  
  labels <- c("Disappearing", "No Change", "Appearing")
  ColScheme <- c(brewer.pal(11,"RdBu")[2], "grey85", brewer.pal(11,"RdBu")[10]); length(ColScheme)
  
  plot(bdy.bc, border="black", lwd=0.4)
  image(X, add=T, xaxt="n", yaxt="n", col=ColScheme, maxpixels= ncell(X))
  legend("topright", legend=labels, fill=ColScheme, bty="n", cex=0.9)
  
  print(spp)
  dev.off()
  
}





