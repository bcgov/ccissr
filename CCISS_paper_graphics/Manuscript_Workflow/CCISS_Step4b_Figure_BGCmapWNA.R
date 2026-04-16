
##======================================================================================
## CCISS Publication Scripts
## Step 4b - Figure - BGC Map for full WNA extent
##======================================================================================


#==================================================
# spatial data
#==================================================

setwd("C:\\Users\\CMAHONY\\OneDrive - Government of BC\\SpatialData\\Boundaries")

# ###country boundaries
# # ORIGIONAL SOURCE: http://www.diva-gis.org/gdata
# countries <- readOGR(dsn="countries", layer='countries')
# countries.NA <- countries[grep("Canada|United States|Mexico", countries$COUNTRY),]
# 
# ####### create a polygon mask for North America.
# my_box = as(extent(-179, -50, -20, 84), "SpatialPolygons")      		# convert extent box to shapefile (rectangle)
# cont.NA <- unionSpatialPolygons(countries.NA, rep(1,length(countries.NA$OBJECTID)))
# cont.NA.g <- gSimplify(cont.NA, tol=0.01, topologyPreserve=TRUE)
# proj4string(my_box) = projection(cont.NA)				# assign spatial projection to extent object
# mask.NA <- gDifference(my_box, cont.NA.g)
# projection(mask.NA)  #verify latlong projection of the study area boundary
# P4S.latlon <- CRS("+proj=longlat +datum=WGS84")

### admin boundaries
bdy.usa1 <- readOGR("USA_adm",'USA_adm1')
bdy.usa <- gSimplify(bdy.usa1, tol=0.01, topologyPreserve=TRUE) #generalize the linework
bdy.can1 <- readOGR("CAN_adm",'CAN_adm1')
bdy.can <- gSimplify(bdy.can1, tol=0.01, topologyPreserve=TRUE) #generalize the linework

# areas <- lapply(bdy.usa@polygons, function(x) sapply(x@Polygons, function(y) y@area))
# bigpolys <- lapply(areas, function(x) which(x > 1000))
# for(i in 1:length(bigpolys)){
#   if(length(bigpolys[[i]]) >= 1 && bigpolys[[i]][1] >= 1){
#     bdy.usa@polygons[[i]]@Polygons <- bdy.usa@polygons[[i]]@Polygons[bigpolys[[i]]]
#     bdy.usa@polygons[[i]]@plotOrder <- 1:length(bdy.usa@polygons[[i]]@Polygons)
#   }
# }
setwd("C:\\GitHub\\2019_CCISS")

#===============================================================================
# Set analysis Parameters
#===============================================================================

source("./_CCISS_Packages.R") ## packages required
source("./_CCISS_Functions.R") ## common functions
source("./_CCISS_Parameters.R") ## settings used through all scripts

# #==================================================
# # BGC projections
# #==================================================

proj.year <- 2055
rcp <- "rcp45"
GCM <- "ensemble"

## parameters
grid <- "WNA2"
BGC.pred.ref <- read.csv(paste("outputs\\BGC.pred", grid, "ref", model, "csv", sep="."), header=F)[,1]
BGC.pred.proj <- read.csv(paste("outputs/BGC.pred", grid, GCM, rcp, proj.year, model, "csv", sep = "."), header=F)[,1]
unique(BGC.pred.ref)

## parameters
grid <- "Salish1"
BGC.pred.ref.inset <- read.csv(paste("outputs\\BGC.pred", grid, "ref", model, "csv", sep="."), header=F)[,1]
BGC.pred.proj.inset <- read.csv(paste("outputs/BGC.pred", grid, GCM, rcp, proj.year, model, "csv", sep = "."), header=F)[,1]
unique(BGC.pred.ref.inset)


##############
# (A) reference BGC map
##############


######################
##reduce subzone-variant to zone

#BGC zone color scheme
BGCcolors$colour <- as.character(BGCcolors$colour)
BGCcolors$colour[match(BGCcolors.BC$zone, BGCcolors$classification)] <- as.character(BGCcolors.BC$HEX)
ColScheme.zone <- factor(BGCcolors$colour, levels=BGCcolors$colour)
zones <- factor(BGCcolors$classification, levels=BGCcolors$classification)

ColScheme.subzone <- BGCcolors.subzone$colour
subzones <- factor(BGCcolors.subzone$classification, levels=BGCcolors.subzone$classification)

zone.pred.ref <- gsub("[:a-z:]","",BGC.pred.ref) 
zone.pred.ref <- gsub("[:1-9:]","",zone.pred.ref) 
zone.pred.ref <- gsub("_.*","",zone.pred.ref)
for(zone in zones){zone.pred.ref[grep(zone,zone.pred.ref)] <- zone}
zone.pred.ref <- factor(zone.pred.ref, levels=zones)

BGC.pred.ref.inset <- factor(as.character(BGC.pred.ref.inset), levels=subzones)

######################
## plot the map

grid <- "WNA2"
grid.dem <- "dem2_WNA"
grid.data <- read.csv(paste("inputs\\", grid, ".csv", sep = ""))
dem <- raster(paste("inputs\\", grid.dem,".tif", sep=""))
land.fine <- which(!is.na(values(dem)))  # raster cells with no dem value
P4S.latlon <- CRS("+proj=longlat +datum=WGS84")

X <- dem
par(mar=c(0,0,0,0))
# plot(X)

values(X) <- NA
values(X)[land.fine] <- zone.pred.ref
# values(X)[1:length(levels(zone.pred.ref))] <- 1:length(levels(zone.pred.ref)) # this is a patch that is necessary to get the color scheme right. 

png(filename=paste("results\\CCISS.manu.BGCmap.ref", model, "png",sep="."), type="cairo", units="in", width=6.5, height=8, pointsize=9, res=600)
# pdf(file=paste("Results\\CCISS.Fig1.BGCmap","pdf",sep="."), width=6.5, height=8, pointsize=9)

par(mar=c(0.1,0.1,0.1,0.1))
# image(hill, xlim=c(-135, -108), ylim=c(39, 60), col=alpha(grey(0:100/100), 1), xaxt="n", yaxt="n", maxpixels= ncell(hill))
image(X, xlim=c(-135, -108), ylim=c(39, 60), xaxt="n", yaxt="n", col=alpha(ColScheme.zone, 1), maxpixels=ncell(X))
# plot(mask.NA, add=T, col="white", border=F)
plot(bdy.usa1, add=T, lwd=0.4)
plot(bdy.can1, add=T, lwd=0.4)

bgcs <- levels(zone.pred.ref)
quantile.lat <- rep(0.5, length(bgcs))
quantile.lat[which(bgcs=="CWH")] <- 0.7
quantile.lat[which(bgcs=="CVG")] <- 0.6
quantile.lat[which(bgcs=="MHRF")] <- 0.6
quantile.lat[which(bgcs=="CWF")] <- 0.48
quantile.lat[which(bgcs=="CMA")] <- 0.3
quantile.lat[which(bgcs=="CRF")] <- 0.55

for(bgc in bgcs){
  i=which(bgcs==bgc)
  temp <- grid.data[which(zone.pred.ref==bgc),]
  select.lat <- temp$lat[max(which(temp$lat >= quantile(temp$lat, quantile.lat[i])))]
  pt <- round(quantile(which(temp$lat==select.lat), 0.5))
  points(temp[pt,4:3], pch=21, bg=alpha(ColScheme.zone[which(bgcs==bgc)], 1), cex=1.5)
  text(temp[pt,4:3]-c(0,0), bgc, pos=4, cex=1, font=2)
  print(paste(which(bgcs==bgc), "-", bgc))
}
box() 
mtext(paste("(a) ", sep=""), side=3, line=-14.5, adj=0.02, cex=1.5, font=2)

grid <- "Salish1"
grid.dem <- "dem1_Salish"
grid.data <- read.csv(paste("inputs\\", grid, ".csv", sep = ""))
dem <- raster(paste("inputs\\", grid.dem,".tif", sep=""))
land.fine <- which(!is.na(values(dem)))  # raster cells with no dem value
length(land.fine)

extent <- c(-124.75, -121,46.5,50)
dem <- crop(dem, extent)
land.fine <- which(!is.na(values(dem)))  # raster cells with no dem value

select.crop <- which(grid.data$lon>extent[1] & grid.data$lon<extent[2] & grid.data$lat>extent[3] & grid.data$lat<extent[4])
grid.data.crop <- grid.data[select.crop, ]

X <- dem
values(X)[land.fine] <- BGC.pred.ref.inset[select.crop]
values(X)[1:length(subzones)] <- 1:length(subzones) # this is a patch that is necessary to get the color scheme right.
xlim=c(extent(X)[1], extent(X)[2])
ylim=c(extent(X)[3], extent(X)[4])
rect(xlim[1],  ylim[1],  xlim[2],  ylim[2],  col=(alpha("white", 0)), lwd=1.5)

par(plt = c(0.01, 0.375, 0.01, 0.45), new = TRUE)
image(X, xaxt="n", yaxt="n", xlim=xlim, ylim=ylim, col=alpha(ColScheme.subzone, 1), maxpixels=ncell(X)) 
# plot(mask.NA, add=T, col="white", border=F)
plot(bdy.usa1, add=T, lwd=0.4)
plot(bdy.can1, add=T, lwd=0.4)
mtext(paste("(b) ", sep=""), side=1, line=-1.5, adj=0.02, cex=1.5, font=2)

bgcs <- subzones
quantile.lat <- rep(0.5, length(bgcs))
quantile.lat[which(bgcs=="CWHvh_WA")] <- 0.35
quantile.lat[which(bgcs=="CWHxm2")] <- 0.89
quantile.lat[which(bgcs=="CDFmm")] <- 0.6
quantile.lat[which(bgcs=="MHmm1")] <- 0.6

for(bgc in bgcs){
  i=which(bgcs==bgc)
  temp <- grid.data.crop[which(BGC.pred.ref.inset[select.crop]==bgc),]
  median.lat <- temp$lat[max(which(temp$lat >= quantile(temp$lat, quantile.lat[i])))]
  pt <- round(quantile(which(temp$lat==median.lat), 0.25))
  points(temp[pt,4:3], pch=21, bg=alpha(ColScheme.subzone[which(bgcs==bgc)], 1), cex=1.5)
  text(temp[pt,4:3]-c(0,0), bgc, pos=4, cex=.8, font=2)
  print(paste(which(bgcs==bgc), "-", bgc))
}

box()
dev.off()


##############
# (B) Projected BGC map
##############

######################
##reduce subzone-variant to zone

zone.pred.proj <- gsub("[:a-z:]","",BGC.pred.proj) 
zone.pred.proj <- gsub("[:1-9:]","",zone.pred.proj) 
zone.pred.proj <- gsub("_.*","",zone.pred.proj)
for(zone in zones){zone.pred.proj[grep(zone,zone.pred.proj)] <- zone}
zone.pred.proj <- factor(zone.pred.proj, levels=zones)

BGC.pred.proj.inset <- factor(as.character(BGC.pred.proj.inset), levels=subzones)

######################
## plot the map

grid <- "WNA2"
grid.dem <- "dem2_WNA"
grid.data <- read.csv(paste("inputs\\", grid, ".csv", sep = ""))
dem <- raster(paste("inputs\\", grid.dem,".tif", sep=""))
land.fine <- which(!is.na(values(dem)))  # raster cells with no dem value
P4S.latlon <- CRS("+proj=longlat +datum=WGS84")

X <- dem
par(mar=c(0,0,0,0))
# plot(X)

values(X) <- NA
values(X)[land.fine] <- zone.pred.proj
# values(X)[1:length(levels(zone.pred.proj))] <- 1:length(levels(zone.pred.proj)) # this is a patch that is necessary to get the color scheme right. 

png(filename=paste("results\\CCISS.manu.BGCmap", GCM, rcp, proj.year, model, "png",sep="."), type="cairo", units="in", width=6.5, height=8, pointsize=9, res=600)
# pdf(file=paste("Results\\CCISS.Fig1.BGCmap","pdf",sep="."), width=6.5, height=8, pointsize=9)

par(mar=c(0.1,0.1,0.1,0.1))
# image(hill, xlim=c(-135, -108), ylim=c(39, 60), col=alpha(grey(0:100/100), 1), xaxt="n", yaxt="n", maxpixels= ncell(hill))
image(X, xlim=c(-135, -108), ylim=c(39, 60), xaxt="n", yaxt="n", col=alpha(ColScheme.zone, 1), maxpixels=ncell(X))
# plot(mask.NA, add=T, col="white", border=F)
plot(bdy.usa1, add=T, lwd=0.4)
plot(bdy.can1, add=T, lwd=0.4)

bgcs <- levels(zone.pred.proj)
quantile.lat <- rep(0.5, length(bgcs))
quantile.lat[which(bgcs=="CWH")] <- 0.7
quantile.lat[which(bgcs=="CVG")] <- 0.6
quantile.lat[which(bgcs=="MHRF")] <- 0.6
quantile.lat[which(bgcs=="CWF")] <- 0.48
quantile.lat[which(bgcs=="CMA")] <- 0.3
quantile.lat[which(bgcs=="CRF")] <- 0.55

for(bgc in bgcs){
  i=which(bgcs==bgc)
  temp <- grid.data[which(zone.pred.proj==bgc),]
  select.lat <- temp$lat[max(which(temp$lat >= quantile(temp$lat, quantile.lat[i])))]
  pt <- round(quantile(which(temp$lat==select.lat), 0.5))
  points(temp[pt,4:3], pch=21, bg=alpha(ColScheme.zone[which(bgcs==bgc)], 1), cex=1.5)
  text(temp[pt,4:3]-c(0,0), bgc, pos=4, cex=1, font=2)
  print(paste(which(bgcs==bgc), "-", bgc))
}
box() 
mtext(paste("(a) ", sep=""), side=3, line=-14.5, adj=0.02, cex=1.5, font=2)

grid <- "Salish1"
grid.dem <- "dem1_Salish"
grid.data <- read.csv(paste("inputs\\", grid, ".csv", sep = ""))
dem <- raster(paste("inputs\\", grid.dem,".tif", sep=""))
land.fine <- which(!is.na(values(dem)))  # raster cells with no dem value
length(land.fine)

extent <- c(-124.75, -121,46.5,50)
dem <- crop(dem, extent)
land.fine <- which(!is.na(values(dem)))  # raster cells with no dem value

select.crop <- which(grid.data$lon>extent[1] & grid.data$lon<extent[2] & grid.data$lat>extent[3] & grid.data$lat<extent[4])
grid.data.crop <- grid.data[select.crop, ]

X <- dem
values(X)[land.fine] <- BGC.pred.proj.inset[select.crop]
values(X)[1:length(subzones)] <- 1:length(subzones) # this is a patch that is necessary to get the color scheme right.
xlim=c(extent(X)[1], extent(X)[2])
ylim=c(extent(X)[3], extent(X)[4])
rect(xlim[1],  ylim[1],  xlim[2],  ylim[2],  col=(alpha("white", 0)), lwd=1.5)

par(plt = c(0.01, 0.375, 0.01, 0.45), new = TRUE)
image(X, xaxt="n", yaxt="n", xlim=xlim, ylim=ylim, col=alpha(ColScheme.subzone, 1), maxpixels=ncell(X)) 
# plot(mask.NA, add=T, col="white", border=F)
plot(bdy.usa1, add=T, lwd=0.4)
plot(bdy.can1, add=T, lwd=0.4)
mtext(paste("(b) ", sep=""), side=1, line=-1.5, adj=0.02, cex=1.5, font=2)

bgcs <- subzones
quantile.lat <- rep(0.5, length(bgcs))
quantile.lat[which(bgcs=="CWHvh_WA")] <- 0.35
quantile.lat[which(bgcs=="CWHxm2")] <- 0.89
quantile.lat[which(bgcs=="CDFmm")] <- 0.6
quantile.lat[which(bgcs=="MHmm1")] <- 0.6

for(bgc in bgcs){
  i=which(bgcs==bgc)
  temp <- grid.data.crop[which(BGC.pred.proj.inset[select.crop]==bgc),]
  median.lat <- temp$lat[max(which(temp$lat >= quantile(temp$lat, quantile.lat[i])))]
  pt <- round(quantile(which(temp$lat==median.lat), 0.25))
  points(temp[pt,4:3], pch=21, bg=alpha(ColScheme.subzone[which(bgcs==bgc)], 1), cex=1.5)
  text(temp[pt,4:3]-c(0,0), bgc, pos=4, cex=.8, font=2)
  print(paste(which(bgcs==bgc), "-", bgc))
}

box()
dev.off()
