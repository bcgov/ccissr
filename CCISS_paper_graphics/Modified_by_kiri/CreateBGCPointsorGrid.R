###This script imports BEC shapefiles and DEMs to create training point data. 
##Process 1 will create a random sample of points per BGC
##Process 2 creates a regular grid of points to be applied in mapping script
###Kiri Daust, August 2018

.libPaths("E:/R packages351")
rm(list=ls())
library(dplyr)
library(rgdal)
library(sp)
library(raster)
library(rgeos)
library(maptools)
library(magrittr)
library(tibble)
library(tidyr)
library(sf)
library(tcltk)
library(foreach)
library(httr)
library(jsonlite)
library(randomForest)
library(data.table)

wd <- tk_choose.dir(); setwd(wd)

dem <- raster("bc25fill") ###Read DEM
bec11 <- st_read(dsn="BGCv11_WithLandcover.gdb",layer="BGCv11_withLandcover") ##read BGC shape file - updated to clipped version
CRS.albers <- CRS ("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs")
allUnits <- unique(as.character(bec11$MAP_LABEL))###What units are in BEC?
allUnits <- allUnits[allUnits !=""]
##set up for parallel processing
require(doParallel)
set.seed(123321)
coreNum <- as.numeric(detectCores()-1)
coreNo <- makeCluster(coreNum)
registerDoParallel(coreNo, cores = coreNum)
clusterEvalQ(coreNo, .libPaths("E:/R packages351"))
#BGC = "CWHwh1"

###randomly select  points within each BEC unit and get elevation data
out <- foreach(BGC = allUnits, .combine = rbind, .packages = c("sf","sp","raster")) %dopar% {
  temp <- bec11$Shape[bec11$MAP_LABEL == BGC] ###Extract polygons for each subzones
  temp <- as(temp, "Spatial") ##conver to sp
  p <- spsample(temp, 1000, type = "random", offset = c(0, 0)) #, iter = 15change for number of points here
  p <- spTransform(p, CRS("+init=epsg:4326")) #to lat lon
  coords <- p@coords
  coords <- as.data.frame(coords)
  
 p2 <- spTransform(p, CRS(proj4string(dem)))
 coords$Elevation <- raster::extract(dem,p2) ##get elevation for each point
  coords$BGC <- BGC
  coords
}

out2 <- out[out$y < 60,]
write.csv(out2,"BECv11_1000Pt.csv", row.names = TRUE) ##rename for version and number of points


####select random points from Alberta
ABdem <- raster("AlbertaDEM.tif")
ABSNR <- st_read(dsn="AlbertaNSR.gdb",layer="AlbertaNR") ##read BGC shape file - updated to clipped version
CRS.albers <- CRS ("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs")
allUnits <- unique(as.character(ABSNR$NSRCODE))###What units are in AB?
allUnits <- allUnits[allUnits !=""]
##set up for parallel processing
require(doParallel)
set.seed(123321)
coreNum <- as.numeric(detectCores()-1)
coreNo <- makeCluster(coreNum)
registerDoParallel(coreNo, cores = coreNum)
clusterEvalQ(coreNo, .libPaths("E:/R packages351"))
#BGC = "NM"
###randomly select 2000 points within each BEC unit and get elevation data
out <- foreach(BGC = allUnits, .combine = rbind, .packages = c("sf","sp","raster")) %dopar% {
  temp <- ABSNR$Shape[ABSNR$NSRCODE == BGC] ###Extract polygons for each subzones
  temp <- as(temp, "Spatial") ##conver to sp
  p <- spsample(temp, 500, type = "random") #, iter = 15change for number of points here
  p <- spTransform(p, CRS("+init=epsg:4326")) #to lat lon
  coords <- p@coords
  coords <- as.data.frame(coords)
  
  p2 <- spTransform(p, CRS(proj4string(ABdem)))
  coords$Elevation <- raster::extract(ABdem,p2) ##get elevation for each point
  coords$BGC <- BGC
  coords
}

out2 <- out[out$y < 60,]
write.csv(out2,"AlbertaSNR_500Pt.csv", row.names = TRUE) ##rename for version and number of points


####4 km grid BC, AB, US########################
CRS.albers <- CRS ("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs")

BC <- readOGR(dsn = "BC_AB_US_Shp", layer = "ProvincialOutline")
AB <- readOGR(dsn = "BC_AB_US_Shp", layer = "AlbertaSubRegions")
US <- readOGR(dsn = "BC_AB_US_Shp", layer = "USA_States")

BC <- spTransform(BC, CRS.albers)
AB <- spTransform(AB, CRS.albers)
US <- spTransform(US, CRS.albers)

BCdem <- raster("bc25fill")

###combine alberta and US DEMs
ABdem <- raster("AlbertaDEM.tif")
USdem1 <- raster("USA_WestDEM.tif")
USdem2 <- raster("USA_WestCoastDEM.tif")
USdem <- raster::merge(USdem1,USdem2)
allDEM <- raster::merge(ABdem,USdem) # not possible to merge rasters due to different resolutions
writeRaster(allDEM, filename = "US_AB_DEM.tif")
allDEM <- raster("US_AB_DEM.tif")

#####Loop through BC, Alberta, and US to create grid
names <- c("BC","AB","USA")
i <-0
grid4k <- foreach(X = c("ProvincialOutline", "AlbertaSubRegions", "USA_States"), .combine = rbind) %do% {
  i <- i+1
  shp <- readOGR(dsn = "BC_AB_US_Shp", layer = X)
  shp <- spTransform(shp, CRS.albers)
  p <- spsample(shp, cellsize = c(1000,1000), type = "regular")# change this to change grid size
  p <- spTransform(p, CRS("+init=epsg:4326")) #to lat long
  coords <- p@coords
  coords <- as.data.frame(coords)
  if(i == 1){
    dem <- BCdem
  }else{
    dem <- allDEM
  }
  p <- spTransform(p, CRS(proj4string(dem)))
  coords$elev <- raster::extract(dem,p)
  colnames(coords) <- c("lon","lat","el")
  coords$Region <- names[i]
  coords
}

rownames(grid4k) <- NULL
setDT(grid4k, keep.rownames = TRUE)[]
grid4k <-grid4k[,c("rn", "Region", "lat", "lon", "el")]
grid4k <- grid4k[grid4k$el >0,]
#grid4US <- grid4US[grid4US$el >0,]
grid4BC <- grid4k[grid4k$Region == "BC",]
grid4USAB <- grid4k[grid4k$Region != "BC",]
grid4USAB <- grid4USAB[grid4USAB$lat >= 37,]
grid4USAB <- grid4USAB[grid4USAB$llon <= -104,]
grid4US <-grid4USAB[grid4USAB$Region != "AB",]
write.csv (grid4k, "BC_AB_US1km.csv", row.names = FALSE)
write.csv(grid4BC, "BC1km.csv", row.names = FALSE)
write.csv(grid4USAB,"US_AB1km.csv", row.names = FALSE)
write.csv (grid4US,  "US1km.csv", row.names = FALSE)


 #############################################################################3
####Assign BGCs and subregions
##################################################################################3
###for BC (input 4 km grid, output is just BC points with BGC assigned)

BCwBGC <- foreach(BGC = allUnits, .combine = rbind) %do%{
  dat <- grid4BC
  pointsOrig <- dat
  coordinates(dat) <- c("lat","lon")
  proj4string(dat) <- CRS("+init=epsg:4326")
  dat <- spTransform(dat, CRS.albers)  # standard albers projection for BC gov't data
  
  tempPoly <- bec11[bec11$MAP_LABEL == BGC,]
  tempPoly <- as(tempPoly, "Spatial") ##conver to sp
  tempPoly <- spTransform(tempPoly, CRS.albers) 
  dat <- over(dat, tempPoly) ###which ones are inside the BGC
  pointsOrig <- pointsOrig[!is.na(dat$BGC_LABEL),] ###Remove points not inside BGC
  if(nrow(pointsOrig) > 0){ ###check that some points fall inside BGC
    pointsOrig$BGC <- BGC
    pointsOrig
  }
}
BCwBGC$ID1 <-row.names(BCwBGC)
colnames(BCwBGC)[5] <- "ID2"
BC4kgrid2  <- BCwBGC [,c("ID1", "ID2","lat", "lon", "el")]
write.csv (BC4kgrid2, "BC2km_BGC.csv")

 ###for AB (input 4 km grid, output is just Alberta points with BGC assigned)
abUnits <- as.character(AB$NSRNAME)

ABwBGC <- foreach(BGC = abUnits, .combine = rbind) %do%{
  dat <- grid4k
  pointsOrig <- dat
  coordinates(dat) <- c("lat","long")
  proj4string(dat) <- CRS("+init=epsg:4326")
  dat <- spTransform(dat, CRS.albers)  # standard albers projection for BC gov't data
  
  tempPoly <- AB[AB$NSRNAME == BGC,]
  tempPoly <- spTransform(tempPoly, CRS.albers) 
  dat <- over(dat, tempPoly) ###which ones are inside the BGC
  pointsOrig <- pointsOrig[!is.na(dat),] ###Remove points not inside BGC
  pointsOrig$BGC <- BGC
  pointsOrig
}

