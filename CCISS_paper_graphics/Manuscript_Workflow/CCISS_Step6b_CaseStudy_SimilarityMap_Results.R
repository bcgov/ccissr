#####Post-processing of climate novelty maps

### NB: THE STRATIFICATION TO GENERATE REPRESENTATIVE LOCATIONS WAS DONE WITH UN-TRANSFORMED DEGREE DAYS AND CMD. THIS WILL NEED TO BE UPDATED FOR THE FINAL ANALYSIS. 
models <-  c("ACCESS1-0","CanESM2","CCSM4","CESM1-CAM5","CNRM-CM5","CSIRO-Mk3-6-0", "GFDL-CM3","GISS-E2R", "GlobalMean","HadGEM2-ES", "INM-CM4", "IPSL-CM5A-MR", "MIROC-ESM", "MIROC5", "MPI-ESM-LR","MRI-CGCM3")    

options(digits=22)  #increase the defualt signif digits to 22 from 7.

library(scales)
library(MASS)   #contains the lda; also eqscplot(x,y) function, which plots axes of equal scale
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
library(rgdal)
library(rgeos)
library(GISTools)
require(adehabitatLT)   #provides calculation of the chi distribution



# #####################
# #### Input data
# #####################
type <- "SimpleMahal"
VarSet <- "Seasonal"
VarCode <- "S"
Grid.fine <- "NAnaec4"
park.analogs <- "SingleParkReverse"
Grid.medium <- "parks.F32"
Grid.analog2 <- ""
Grid.analog3 <- ""
Grid.map <- "dem4"
proj.year <- 2085
model.fine <- if(Grid.fine%in%c("parks.F16", "parks.F32", "parks.F64")) "15GCM-Ensemble" else "GlobalMean"
RCP <- "RCP45"
scenario <- "6variable"


## DEM and land definition
setwd("D:\\Backup_Jun12\\Masters\\Research\\Data\\BCparks\\InputData")
dem <- raster(paste(Grid.map,".tif", sep=""))
land.medium <- read.csv(paste("land.",Grid.medium,".csv", sep=""))[,1]  
land.fine <- read.csv(paste("land.",Grid.fine,".csv", sep=""))[,1]  

#read in identities and distances to 4 nearest CRU neighbours for each reference grid cell. 
A.grid.4nn <- read.csv(paste("grid4nn_",Grid.fine,".csv",sep=""), strip.white = TRUE, na.strings = c("NA","",-9999) )
A.grid.dist <- read.csv(paste("griddist_",Grid.fine,".csv",sep=""), strip.white = TRUE, na.strings = c("NA","",-9999) )

trunc.SD2s <- 0.1  
trunc.rule2 <- paste(trunc.SD2s,"SD2",sep="")

# setwd("D:\\Backup_Jun12\\Masters\\Research\\Data\\BCparks\\InputData")
# nonCNA.medium <- read.csv(paste("nonCNA.",Grid.medium,".csv", sep=""))[,1]
# subsample.medium <- read.csv(paste("subsample.",Grid.medium,".csv", sep=""))[,1]
# pts <- read.csv(paste(Grid.medium,".csv", sep=""))[-nonCNA.medium,]
# parks <- read.csv(paste("BCparksSpatial.",Grid.medium,".csv", sep=""))[-nonCNA.medium,]
# park <- pts[grep(selected.park,parks$PROT_NAME),]

################
### Other spatial data
################
setwd("D:\\Backup_Jun12\\Masters\\Research\\SpatialData\\Boundaries")

###country boundaries
# ORIGIONAL SOURCE: http://www.diva-gis.org/gdata
countries <- readOGR(dsn="countries", layer='countries')
countries.NA <- countries[grep("Canada|United States|Mexico", countries$COUNTRY),]

# plot(countries[grep("Canada|United States|Mexico", countries$COUNTRY),])
# plot(countries.NA,xlim=c(-125.75,-125.5),ylim=c(49,49.5))

####### create a polygon mask for North America.
my_box = as(extent(-179, -50, -20, 84), "SpatialPolygons")					# convert extent box to shapefile (rectangle)
proj4string(my_box) = projection(countries.NA)				# assign spatial projection to extent object
mask.NA <- gDifference(my_box, countries.NA)
# plot(mask.NA)
projection(mask.NA)  #verify latlong projection of the study area boundary
P4S.NAEC <- CRS(projection(dem))   #establish the projection of the dem
mask.NAEC <- spTransform(mask.NA, P4S.NAEC) #reproject the countries polygons
#
### admin boundaries
bdy.can1 <- readOGR("CAN_adm",'CAN_adm1')
bdy.usa1 <- readOGR("USA_adm",'USA_adm1')
# bdy.mex1 <- readOGR("MEX_adm",'MEX_adm1')
bdy.can1 <- spTransform(bdy.can1, P4S.NAEC) #reproject the countries polygons
bdy.usa1 <- spTransform(bdy.usa1, P4S.NAEC) #reproject the countries polygons
# bdy.mex1 <- spTransform(bdy.mex1, P4S.NAEC) #reproject the countries polygons

## WWF terrestrial ecoregions of the world (for descriptive purposes)
teow <- readShapePoly("D:\\Backup_Jun12\\Masters\\Research\\SpatialData\\WWF\\TEOW\\wwf_terr_ecos.shp")
P4S.latlon <- CRS("+proj=longlat +datum=WGS84")
projection(teow) <- P4S.latlon
teow <- spTransform(teow, P4S.NAEC)

### BC Parks
P4S.AEA <- CRS ("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs") # standard albers projection for BC gov't data
parks <- readShapePoly("D:\\Backup_Jun12\\Masters\\Research\\SpatialData\\BCparks\\combined\\ParksConservancies.shp")
projection(parks) <- P4S.AEA # tell R the projection of the shapefile (this is known via the shapefile metadata)
# P4S.latlon <- CRS(projection(dem))
# parks.latlon <- spTransform(parks, P4S.latlon) # reproject to lat-long
parks.naec <- spTransform(parks, P4S.NAEC) # reproject to lat-long

# 
# setwd(paste("D:\\Backup_Jun12\\Masters\\Research\\Data\\BCparks\\Results\\", nickname, sep=""))
# png(filename=paste("SimilarityMap",nickname,Grid.fine,Grid.medium,model.fine,RCP,proj.year,trunc.rule2,type,scenario,".png",sep="_"), type="cairo", units="in", width=7.5, height=5, pointsize=10, res=800)
# 
# par(mar=c(0.1,0.1,0.1,0.1), mfrow=c(1,2))
# # par(mgp=c(2.5,0.75,0))
# xl <- -3000000; yb <- -1000000; xr <- -2800000; yt <- 500000
# breakseq <- c(0,2,4)
# breakpoints <- c(seq(breakseq[1], breakseq[3],0.001),99); length(breakpoints)
# ColScheme <- matlab.like2(length(breakpoints)-1)
# ColScheme <- colorRampPalette(c("#0000CD", blue2red(length(breakpoints)-1), "#CD0000"))(length(breakpoints)-1)
# ColScheme <- colorRampPalette(c("light gray", "black"))(length(breakpoints)-1)
# ColScheme <- colorRampPalette(c("yellow", "gray70", "gray30", "black"))(length(breakpoints)-1)
# ColScheme <- colorRampPalette(c("yellow", "gray10", "gray40", "gray70"))(length(breakpoints)-1)
# X <- dem
# 
# for(direction in c("forward", "backward")){
# 
# #import the novelty distance files for the ensemble
#   setwd(paste("D:\\Backup_Jun12\\Masters\\Research\\Data\\BCparks\\OutputData\\", nickname, sep=""))
#   IDWMean <- function(x){IDWAvg <- sum(x[1:(length(x)/2)]/x[(length(x)/2+1):length(x)])/sum(1/x[(length(x)/2+1):length(x)]); return(IDWAvg)}   #this is an awkward solution to the problem of needing a one-argument function for use by "apply". the idea is to cbind the NN.dist.proj and grid.dist data frames into one argument for the function to operate on as a single vector.
# filelist <- list.files(pattern = paste(if(direction == "forward") "NNdistref" else "NNdistproj",Grid.fine,park.analogs , Grid.medium,Grid.analog2,Grid.analog3,model.fine,RCP,proj.year,trunc.rule2,type,nickname,scenario,".csv", sep="_"), recursive =F )
# NN.dist.proj <- read.csv(filelist[grep(model.fine,filelist)])
# NN.dist.proj <- apply(cbind(NN.dist.proj[1:dim(A.grid.dist)[1],],A.grid.dist),1,IDWMean) #apply inverse distance-weighted averaging function
# dimensions <- 6 # i did it this way for simplicity, but it assumes that all measurements were made with 6 pcs. 
# NN.chi.proj <- qchi(pchisq(NN.dist.proj^2,dimensions),1)
# NN.chi.proj[!is.finite(NN.chi.proj)] <- max(NN.chi.proj[is.finite(NN.chi.proj)])
# setwd(paste("D:\\Backup_Jun12\\Masters\\Research\\Data\\BCparks\\Results\\", nickname, sep=""))
# 
# ## distance
# values(X) <- NA
# values(X)[land.fine] <- NN.chi.proj
# plot(mask.NAEC, xaxt="n", yaxt="n", xlim=c(-3500000, -300000), ylim=c(-1000000, 3200000), col="white", border=F) 
# plot(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks=breakpoints, legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 
# plot(bdy.can1, border=alpha("white", 0.65), lwd=0.4, add=T)
# plot(bdy.usa1, border=alpha("white", 0.65), lwd=0.4, add=T)
# plot(mask.NAEC, add=T, col="white", border=F)
# # plot(dem4points, add=T, pch=16, cex=0.001)
# rect(xl-15000,  yb-3000,  xr+14000,  yt+3000,  col=(alpha("white", 1)), border=NA)
# rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme, border=NA)
# rect(xl,  yb,  xr,  yt,  col=NA)
# #   stats <- c(round(cellStats(X,"min", na.rm=T),2), round(mean(c(cellStats(X,"max"), cellStats(X,"min"))),2), round(cellStats(X,"max"),2))
# text(rep(xr,3),c(yb,mean(c(yt,yb)),yt-100000),sapply(c(bquote(.(breakseq[1])*sigma), bquote(.(breakseq[2])*sigma),bquote(.(breakseq[3])*sigma)),as.expression),pos=4,cex=0.9,font=1)  
# text(xl-100000,mean(c(yb,yt)), "Climatic Similarity", font=1, cex=0.9, srt=90) 
# # rect(xl,  yt+100000,  xr,  yt+300000,  col=ColScheme[length(ColScheme)])
# # text(xr,  yt+200000,  bquote(">"*.(breakseq[3])),pos=4,cex=1,font=1)  
# 
# # legend("topleft", legend=description, bty="n", cex=0.9, inset=-0.02)
# # plot(mask.NAEC, add=T, col="white", border=F)
# if(direction=="backward"){
#   text(xr, 1800000, paste("Where the climates\nof ", title, "\nare going", sep=""), font=2, cex=0.9)
#   text(xr, 1100000, paste("Projected future analogs\nfor the historical climates\nof ", title, "\n(",if(RCP=="RCP45") "RCP4.5, " else "RCP8.5, ", substr(proj.year, 1,3),"0s)",  sep=""), font=1, cex=0.8)
# } else {
#   text(xr, 1800000, paste("Where the climates\nof ", title, "\nare coming from", sep=""), font=2, cex=0.9)
#   text(xr, 1100000, paste("Historical analogs for the\nprojected future climates\nof ", title, "\n(",if(RCP=="RCP45") "RCP4.5, " else "RCP8.5, ", substr(proj.year, 1,3),"0s)",  sep=""), font=1, cex=0.8)
# }
# centroid <- getSpPPolygonsLabptSlots(parks.naec[grep(selected.park,parks.naec$PROT_NAME)[1],])
# points(centroid, pch = 1, col = "black", cex=1.9)
# plot(parks.naec[grep(selected.park,parks.naec$PROT_NAME),], border="black", add=T, lwd=0.4)
# text(centroid, title, pos=4, font=2, col="black", cex=0.7)
# box(col="black", lwd=1.5)
# }
# dev.off()

#####################
## Results loop
#####################

setwd("D:\\Backup_Jun12\\Masters\\Research\\Data\\BCparks\\InputData")
complex.lookup <- read.csv(paste("ComplexLookup.csv", sep=""), stringsAsFactors = F, strip.white = T)
nicknames <- unique(complex.lookup$nickname)

nickname=nicknames[8]
for(nickname in nicknames){
  complex.parks <- complex.lookup$PROT_NAME[which(complex.lookup$nickname==nickname)]


#import the novelty distance files for the ensemble
setwd(paste("D:\\Backup_Jun12\\Masters\\Research\\Data\\BCparks\\OutputData\\", nickname, sep=""))
IDWMean <- function(x){IDWAvg <- sum(x[1:(length(x)/2)]/x[(length(x)/2+1):length(x)])/sum(1/x[(length(x)/2+1):length(x)]); return(IDWAvg)}   #this is an awkward solution to the problem of needing a one-argument function for use by "apply". the idea is to cbind the NN.dist.proj and grid.dist data frames into one argument for the function to operate on as a single vector.
filelist <- list.files(pattern = paste("NNdistproj",Grid.fine,park.analogs , Grid.medium,Grid.analog2,Grid.analog3,model.fine,RCP,proj.year,trunc.rule2,type,nickname,scenario,".csv", sep="_"), recursive =F )
NN.dist.proj <- read.csv(filelist[grep(model.fine,filelist)])
NN.dist.proj <- apply(cbind(NN.dist.proj[1:dim(A.grid.dist)[1],],A.grid.dist),1,IDWMean) #apply inverse distance-weighted averaging function
filelist <- list.files(pattern = paste("NNdistref",Grid.fine,park.analogs , Grid.medium,Grid.analog2,Grid.analog3,model.fine,RCP,proj.year,trunc.rule2,type,nickname,scenario,".csv", sep="_"), recursive =F )
NN.dist.ref <- read.csv(filelist[grep(model.fine,filelist)])
NN.dist.ref <- apply(cbind(NN.dist.ref[1:dim(A.grid.dist)[1],],A.grid.dist),1,IDWMean) #apply inverse distance-weighted averaging function
filelist <- list.files(pattern = paste("NNdisthis",Grid.fine,park.analogs , Grid.medium,Grid.analog2,Grid.analog3,model.fine,RCP,proj.year,trunc.rule2,type,nickname,scenario,".csv", sep="_"), recursive =F )
NN.dist.his <- read.csv(filelist[grep(model.fine,filelist)])
NN.dist.his <- apply(cbind(NN.dist.his[1:dim(A.grid.dist)[1],],A.grid.dist),1,IDWMean) #apply inverse distance-weighted averaging function



####################
## individual maps

# # Western North America map settings
# xl <- -3000000; yb <- -1000000; xr <- -2800000; yt <- 500000 #Zoomed out
# x1 <- 15000; y1 <- 3000; x2 <- 14000; y2 <- 3000; x3 <- 100000; y3 <- 100000; y4 <- 1800000; y5 <- 1100000  #Zoomed out
# xlim=c(-3500000, -300000); ylim=c(-1000000, 3200000) #Zoomed out

# NorthWestern North America map settings
xl <- -2700000; yb <- -100000; xr <- -2600000; yt <- 1000000
x1 <- 15000; y1 <- 3000; x2 <- 14000; y2 <- 3000; x3 <- 80000; y3 <- 0; y4 <- 1900000; y5 <- 1400000  #Zoomed out
xlim=c(-3100000, -800000); ylim=c(-100000, 2700000) #Coastal

# # Coastal map settings
# xl <- -2700000; yb <- 100000; xr <- -2600000; yt <- 1000000
# x1 <- 15000; y1 <- 3000; x2 <- 14000; y2 <- 3000; x3 <- 60000; y3 <- 100000; y4 <- 1900000; y5 <- 1400000  #Zoomed out
# xlim=c(-3300000, -1200000); ylim=c(100000, 2700000) #Coastal

breakseq <- c(0,2.5,5)
breakpoints <- c(seq(breakseq[1], breakseq[3],0.05),99); length(breakpoints)
# ColScheme <- matlab.like2(length(breakpoints)-1)
# ColScheme <- colorRampPalette(c("#0000CD", blue2red(length(breakpoints)-1), "#CD0000"))(length(breakpoints)-1)
# ColScheme <- colorRampPalette(c("light gray", "black"))(length(breakpoints)-1)
ColScheme <- colorRampPalette(c("yellow", "yellow", "gray20", "gray35", "gray50", "gray75", "gray80"))(length(breakpoints)-1)
# ColScheme <- colorRampPalette(c("blue", "yellow", "gray10", "gray40", "gray70"))(length(breakpoints)-1)
# ColScheme <- colorRampPalette(c("yellow", "yellow", "blue", "gray30", "gray50", "gray70", "gray90"))(length(breakpoints)-1)
X <- dem

slots <- getSpPPolygonsLabptSlots(parks.naec[which(parks.naec$PROT_NAME%in%complex.parks)[1],])
centroid <- slots[which(slots[,1]==max(slots[,1])),]

direction="historical"
for(direction in c("forward", "historical", "backward")){
  

# x11(width=6.5, height=8.5, pointsize=16)
setwd(paste("D:\\Backup_Jun12\\Masters\\Research\\Data\\BCparks\\Results\\", nickname, sep=""))
  png(filename=paste("SimilarityMap_Mdist",direction,nickname,Grid.fine,Grid.medium,model.fine,RCP,proj.year,trunc.rule2,type,scenario,".png",sep="_"), type="cairo", units="in", width=6.5, height=8.5, pointsize=16, res=600)
  
  par(mar=c(0.1,0.1,0.1,0.1), mfrow=c(1,1))
  
  ## distance
  values(X) <- NA
  values(X)[land.fine] <- if(direction == "forward") NN.dist.proj else if(direction == "backward") NN.dist.ref else NN.dist.his
  plot(mask.NAEC, xaxt="n", yaxt="n", xlim=xlim, ylim=ylim, col="white", border=F) 
  plot(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks=breakpoints, legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 
  plot(bdy.can1, border=alpha("white", 0.65), lwd=0.5, add=T)
  plot(bdy.usa1, border=alpha("white", 0.65), lwd=0.5, add=T)
  plot(mask.NAEC, add=T, col="white", border=F)
  plot(parks.naec, border=alpha("dodgerblue", 0.7), lwd=0.5, add=T)
  rect(xl-x1,  yb-y1,  xr+x2,  yt+y2,  col=(alpha("white", 1)), border=NA)
  rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme, border=NA)
  rect(xl,  yb,  xr,  yt,  col=NA)
  #   stats <- c(round(cellStats(X,"min", na.rm=T),2), round(mean(c(cellStats(X,"max"), cellStats(X,"min"))),2), round(cellStats(X,"max"),2))
  text(rep(xr,3),c(yb,mean(c(yt,yb)),yt-y3),sapply(c(bquote(.(breakseq[1])), bquote(.(breakseq[2])),bquote(.(breakseq[3]))),as.expression),pos=4,cex=0.9,font=1)  
  text(xl-x3,mean(c(yb,yt)), "Climatic distance", font=1, cex=0.9, srt=90) 
  # rect(xl,  yt+100000,  xr,  yt+300000,  col=ColScheme[length(ColScheme)])
  # text(xr,  yt+200000,  bquote(">"*.(breakseq[3])),pos=4,cex=1,font=1)  
  
  # legend("topleft", legend=description, bty="n", cex=0.9, inset=-0.02)
  # plot(mask.NAEC, add=T, col="white", border=F)
  title <- paste(unique(complex.lookup$Complex[which(complex.lookup$nickname==nickname)]), "complex")
  if(direction=="forward"){
    text(xr, y4, paste("Where the historical\nclimates of the\n", title, "\nare going", sep=""), font=2, cex=0.9)
    text(xr, y5, paste("Projected future analogs\nfor the historical climates\nof the ", title, "\n(",if(RCP=="RCP45") "RCP4.5, " else "RCP8.5, ", substr(proj.year, 1,3),"0s)",  sep=""), font=1, cex=0.8)
  } 
  if(direction=="historical"){
    text(xr, y4, paste("Historical climates of the\n", title, sep=""), font=2, cex=0.9)
    text(xr, y5, paste("Historical climatic similarity\nto the historical climates of\nthe ", title,  sep=""), font=1, cex=0.8)
  }
  if(direction=="backward"){
    text(xr, y4, paste("Where the future\nclimates of the\n", title, "\nare coming from", sep=""), font=2, cex=0.9)
    text(xr, y5, paste("Historical analogs for the\nprojected future climates\nof the ", title, "\n(",if(RCP=="RCP45") "RCP4.5, " else "RCP8.5, ", substr(proj.year, 1,3),"0s)",  sep=""), font=1, cex=0.8)
  }
  
  title2 <- paste(unique(complex.lookup$Complex[which(complex.lookup$nickname==nickname)]), "\ncomplex")
  plot(parks.naec[which(parks.naec$PROT_NAME%in%complex.parks),], border="black", add=T, lwd=1)
  text(centroid[1], centroid[2], title2, pos=4, font=2, col="black", cex=0.7)
  box(col="black", lwd=1.5)
dev.off()

print(direction)
}


setwd(paste("D:\\Backup_Jun12\\Masters\\Research\\Data\\BCparks\\Results\\", nickname, sep=""))
png(filename=paste("SimilarityMap_Mdist",nickname,Grid.fine,Grid.medium,model.fine,RCP,proj.year,trunc.rule2,type,scenario,".png",sep="_"), type="cairo", units="in", width=7.5, height=5, pointsize=10, res=800)

par(mar=c(0.1,0.1,0.1,0.1), mfrow=c(1,2))

for(direction in c("forward", "backward")){

  ## distance
  values(X) <- NA
  values(X)[land.fine] <- if(direction == "forward") NN.dist.proj else if(direction == "backward") NN.dist.ref else NN.dist.his
  plot(mask.NAEC, xaxt="n", yaxt="n", xlim=xlim, ylim=ylim, col="white", border=F) 
  plot(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks=breakpoints, legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 
  plot(bdy.can1, border=alpha("white", 0.65), lwd=0.25, add=T)
  plot(bdy.usa1, border=alpha("white", 0.65), lwd=0.25, add=T)
  plot(mask.NAEC, add=T, col="white", border=F)
  plot(parks.naec, border=alpha("dodgerblue", 0.7), lwd=0.5, add=T)
  rect(xl-x1,  yb-y1,  xr+x2,  yt+y2,  col=(alpha("white", 1)), border=NA)
  rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme, border=NA)
  rect(xl,  yb,  xr,  yt,  col=NA)
  #   stats <- c(round(cellStats(X,"min", na.rm=T),2), round(mean(c(cellStats(X,"max"), cellStats(X,"min"))),2), round(cellStats(X,"max"),2))
  text(rep(xr,3),c(yb,mean(c(yt,yb)),yt-y3),sapply(c(bquote(.(breakseq[1])), bquote(.(breakseq[2])),bquote(.(breakseq[3]))),as.expression),pos=4,cex=0.9,font=1)  
  text(xl-x3,mean(c(yb,yt)), "Climatic distance", font=1, cex=0.9, srt=90) 
  # rect(xl,  yt+100000,  xr,  yt+300000,  col=ColScheme[length(ColScheme)])
  # text(xr,  yt+200000,  bquote(">"*.(breakseq[3])),pos=4,cex=1,font=1)  
  
  # legend("topleft", legend=description, bty="n", cex=0.9, inset=-0.02)
  # plot(mask.NAEC, add=T, col="white", border=F)
  title <- paste(unique(complex.lookup$Complex[which(complex.lookup$nickname==nickname)]), "complex")
  if(direction=="forward"){
    text(xr, y4, paste("Where the historical\nclimates of the\n", title, "\nare going", sep=""), font=2, cex=0.9)
    text(xr, y5, paste("Projected future analogs\nfor the historical climates\nof the ", title, "\n(",if(RCP=="RCP45") "RCP4.5, " else "RCP8.5, ", substr(proj.year, 1,3),"0s)",  sep=""), font=1, cex=0.8)
  } 
  if(direction=="backward"){
    text(xr, y4, paste("Where the future\nclimates of the\n", title, "\nare coming from", sep=""), font=2, cex=0.9)
    text(xr, y5, paste("Historical analogs for the\nprojected future climates\nof the ", title, "\n(",if(RCP=="RCP45") "RCP4.5, " else "RCP8.5, ", substr(proj.year, 1,3),"0s)",  sep=""), font=1, cex=0.8)
  }

  title2 <- paste(unique(complex.lookup$Complex[which(complex.lookup$nickname==nickname)]), "\ncomplex")
  plot(parks.naec[which(parks.naec$PROT_NAME%in%complex.parks),], border="black", add=T, lwd=1)
  text(centroid[1], centroid[2], title2, pos=4, font=2, col="black", cex=0.7)
  box(col="black", lwd=1.5)

  print(direction)
}
dev.off()


#simplified 3-panel map
setwd(paste("D:\\Backup_Jun12\\Masters\\Research\\Data\\BCparks\\Results\\", nickname, sep=""))
png(filename=paste("SimilarityMap_Mdist_3panel",nickname,Grid.fine,Grid.medium,model.fine,RCP,proj.year,trunc.rule2,type,scenario,".png",sep="_"), type="cairo", units="in", width=7.5*3/2, height=5, pointsize=18, res=800)
# pdf(filename=paste("SimilarityMap_Mdist_3panel",nickname,Grid.fine,Grid.medium,model.fine,RCP,proj.year,trunc.rule2,type,scenario,".png",sep="_"), width=7.5*3/2, height=5, pointsize=14)

par(mar=c(0.1,0.1,0.1,0.1), mfrow=c(1,3))

# NorthWestern North America map settings
xl <- -2700000; yb <- -100000; xr <- -2600000; yt <- 1000000
x1 <- 15000; y1 <- 3000; x2 <- 14000; y2 <- 3000; x3 <- 80000; y3 <- 0; y4 <- 1900000; y5 <- 1400000  #Zoomed out
xlim=c(-3100000, -800000); ylim=c(-100000, 2700000) #Coastal

for(direction in c("historical", "backward", "forward")){
  
  ## distance
  values(X) <- NA
  values(X)[land.fine] <- if(direction == "forward") NN.dist.proj else if(direction == "backward") NN.dist.ref else NN.dist.his
  plot(mask.NAEC, xaxt="n", yaxt="n", xlim=xlim, ylim=ylim, col="white", border=F) 
  plot(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks=breakpoints, legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 
  plot(bdy.can1, border=alpha("white", 0.65), lwd=0.25, add=T)
  plot(bdy.usa1, border=alpha("white", 0.65), lwd=0.25, add=T)
  plot(mask.NAEC, add=T, col="white", border=F)
  plot(parks.naec, border=alpha("dodgerblue", 0.7), lwd=0.5, add=T)
  if(direction=="historical"){
    rect(xl-x1,  yb-y1,  xr+x2,  yt+y2,  col=(alpha("white", 1)), border=NA)
  rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme, border=NA)
  rect(xl,  yb,  xr,  yt,  col=NA)
  #   stats <- c(round(cellStats(X,"min", na.rm=T),2), round(mean(c(cellStats(X,"max"), cellStats(X,"min"))),2), round(cellStats(X,"max"),2))
  text(rep(xr,3),c(yb,mean(c(yt,yb)),yt-y3),c("High similarity", "Low similarity", "No similarity"),pos=4,cex=0.9,font=1)
  # text(rep(xr,3),c(yb,mean(c(yt,yb)),yt-y3),c("Good match", "Poor match", "No match"),pos=4,cex=0.9,font=1)
  # text(rep(xr,3),c(yb,mean(c(yt,yb)),yt-y3),c("Good climate analog", "Poor climate analog", "No climate analog"),pos=4,cex=0.9,font=1)
  text(xl-x3,mean(c(yb,yt)), "Climatic dissimilarity", font=2, cex=1, srt=90) 
  # rect(xl,  yt+100000,  xr,  yt+300000,  col=ColScheme[length(ColScheme)])
  # text(xr,  yt+200000,  bquote(">"*.(breakseq[3])),pos=4,cex=1,font=1)  
  }
  
  # legend("topleft", legend=description, bty="n", cex=0.9, inset=-0.02)
  # plot(mask.NAEC, add=T, col="white", border=F)
  title <- paste(unique(complex.lookup$Complex[which(complex.lookup$nickname==nickname)]), "complex")
  if(direction=="historical"){
    text(xr, y4, paste("Historical climates\nof the\n", title, sep=""), font=2, cex=0.9)
    text(xr, y5, paste("Historical climatic similarity\nto the historical climates of\nthe ", title,  sep=""), font=1, cex=0.8)
  }
  if(direction=="forward"){
    text(xr, y4, paste("Where the historical\nclimates of the\n", title, "\nare going", sep=""), font=2, cex=0.9)
    text(xr, y5, paste("Projected future analogs\nfor the historical climates\nof the ", title, "\n(",if(RCP=="RCP45") "RCP4.5, " else "RCP8.5, ", substr(proj.year, 1,3),"0s)",  sep=""), font=1, cex=0.8)
  } 
  if(direction=="backward"){
    text(xr, y4, paste("Where the future\nclimates of the\n", title, "\nare coming from", sep=""), font=2, cex=0.9)
    text(xr, y5, paste("Historical analogs for the\nprojected future climates\nof the ", title, "\n(",if(RCP=="RCP45") "RCP4.5, " else "RCP8.5, ", substr(proj.year, 1,3),"0s)",  sep=""), font=1, cex=0.8)
  }
  
  title2 <- paste(unique(complex.lookup$Complex[which(complex.lookup$nickname==nickname)]), "\ncomplex")
  plot(parks.naec[which(parks.naec$PROT_NAME%in%complex.parks),], border="black", add=T, lwd=1)
  text(centroid[1], centroid[2], title2, pos=4, font=2, col="black", cex=0.7)
  box(col="black", lwd=1.5)
  
  print(direction)
}
dev.off()

print(nickname)
}