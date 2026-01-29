
##======================================================================================
## CCISS Publication Scripts
## Step 4b - Figure - BGC Map for full WNA extent
##======================================================================================

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


#==================================================
# spatial data
#==================================================


setwd("C:\\Colin\\SpatialData\\Boundaries")

###country boundaries
# ORIGIONAL SOURCE: http://www.diva-gis.org/gdata
countries <- readOGR(dsn="countries", layer='countries')
countries.NA <- countries[grep("Canada|United States|Mexico", countries$COUNTRY),]

####### create a polygon mask for North America.
my_box = as(extent(-179, -50, -20, 84), "SpatialPolygons")      		# convert extent box to shapefile (rectangle)
cont.NA <- unionSpatialPolygons(countries.NA, rep(1,length(countries.NA$OBJECTID)))
cont.NA.g <- gSimplify(cont.NA, tol=0.01, topologyPreserve=TRUE)
proj4string(my_box) = projection(cont.NA)				# assign spatial projection to extent object
mask.NA <- gDifference(my_box, cont.NA.g)
projection(mask.NA)  #verify latlong projection of the study area boundary
P4S.latlon <- CRS("+proj=longlat +datum=WGS84")

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
setwd("C:\\Colin\\Projects\\2019_CCISS")


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
spps.lookup <- read.csv("InputData\\Tree speciesand codes_2.0.csv")
edatope.name <- c("Subxeric-poor", "Mesic-medium", "Hygric-rich")
proj.year.name=c("2020s", "2050s", "2080s")
rcp.name=c("RCP4.5", "RCP8.5")
BGCcolors <- read.csv("C:\\Colin\\Projects\\2019_CCISS\\InputData\\BGCzone_Colorscheme.csv")


#===============================================================================
# BGC change for each model prediction
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


# #==================================================
# # BGC projections
# #==================================================
setwd("C:\\Colin\\Projects\\2019_CCISS")

## parameters
grid <- "WNA2"
grid.dem <- "dem2_WNA"

###Load random forest model
fname="Rcode\\BGCv11_AB_USA_16VAR_SubZone_RFmodel.Rdata"
#fname = (file.choose())
load(fname)


fplot=paste("InputData\\", grid, "_Normal_1961_1990MSY.csv", sep="")

Columns <- unique(c("PPT05", "PPT06", "PPT07", "PPT08", "PPT09", "PPT_at", "PPT_wt", "CMD07", "CMD", rownames(importance(BGCmodel))[-which(rownames(importance(BGCmodel))%in%c("PPT_MJ", "PPT_JAS", "PPT.dormant", "CMD.def", "CMDMax", "CMD.total"))]))

Y0 <- fread(fplot, select = Columns, stringsAsFactors = FALSE, data.table = FALSE) #fread is faster than read.csv

#####generate some additional variables
Y0$PPT_MJ <- Y0$PPT05 + Y0$PPT06 # MaY/June precip
Y0$PPT_JAS <- Y0$PPT07 + Y0$PPT08 + Y0$PPT09 # July/Aug/Sept precip
Y0$PPT.dormant <- Y0$PPT_at + Y0$PPT_wt # for calculating spring deficit
Y0$CMD.def <- 500 - (Y0$PPT.dormant)# start of growing season deficit original value was 400 but 500 seems better
Y0$CMD.def [Y0$CMD.def < 0] <- 0 #negative values set to zero = no deficit
Y0$CMDMax <- Y0$CMD07
Y0$CMD.total <- Y0$CMD.def + Y0$CMD

##Predict reference subzones###### (need to break it up because it is big)
BGC.pred.ref <- vector()
factor <- 9
for(i in 1:factor){
  temp <- Y0[(dim(Y0)[1]/factor*(i-1)+1):(dim(Y0)[1]/factor*i),]
  # str(temp)
  BGC.pred.ref <- c(BGC.pred.ref, as.character(predict(BGCmodel, temp)))
  print(i)
}
BGC.pred.ref
write.csv(BGC.pred.ref, paste("OutputData\\BGC.pred.ref", grid, "csv", sep="."), row.names = F)

BGC.pred.ref <- read.csv(paste("OutputData\\BGC.pred.ref", grid, "csv", sep="."))[,1]
unique(BGC.pred.ref)

# #==================================================
# # BGC projections
# #==================================================
setwd("C:\\Colin\\Projects\\2019_CCISS")

## parameters
grid <- "Salish1"
grid.dem <- "dem1_Salish"

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

##Predict reference subzones###### (need to break it up because it is big)
BGC.pred.ref.inset <- as.character(predict(BGCmodel, Y0))
# BGC.pred.ref
write.csv(BGC.pred.ref.inset, paste("OutputData\\BGC.pred.ref", grid, "csv", sep="."), row.names = F)

BGC.pred.ref.inset <- read.csv(paste("OutputData\\BGC.pred.ref", grid, "csv", sep="."))[,1]
unique(BGC.pred.ref.inset)

#===============================================================================
# Projected BGC zones
#===============================================================================

## parameters
grid <- "WNA2"
grid.dem <- "dem2_WNA"
rcp="rcp45"
proj.year=2055

# fplot=paste("InputData\\WNA2_15GCM-Ensemble_rcp45_2055MSY.csv", sep="")
# 
# Y0 <- fread(fplot, select = Columns, stringsAsFactors = FALSE, data.table = FALSE) #fread is faster than read.csv
# 
# Y0 <- Y0[!is.na(Y0[,2]),]
# 
# #####generate some additional variables
# Y0$PPT_MJ <- Y0$PPT05 + Y0$PPT06 # MaY/June precip
# Y0$PPT_JAS <- Y0$PPT07 + Y0$PPT08 + Y0$PPT09 # July/Aug/Sept precip
# Y0$PPT.dormant <- Y0$PPT_at + Y0$PPT_wt # for calculating spring deficit
# Y0$CMD.def <- 500 - (Y0$PPT.dormant)# start of growing season deficit original value was 400 but 500 seems better
# Y0$CMD.def [Y0$CMD.def < 0] <- 0 #negative values set to zero = no deficit
# Y0$CMDMax <- Y0$CMD07
# Y0$CMD.total <- Y0$CMD.def + Y0$CMD
# 
# ##Predict future subzones######
# BGC.pred <- vector()
# factor <- 9
# for(i in 1:factor){
#   temp <- Y0[(dim(Y0)[1]/factor*(i-1)+1):(dim(Y0)[1]/factor*i),]
#   # str(temp)
#   BGC.pred <- c(BGC.pred, as.character(predict(BGCmodel, temp)))
#   print(i)
# }
# write.csv(BGC.pred, paste("OutputData\\BGC.pred", grid, rcp, proj.year, "csv", sep="."), row.names = F)

BGC.pred <- read.csv(paste("OutputData\\BGC.pred", grid, rcp, proj.year, "csv", sep="."))[,1]
unique(BGC.pred)


#===============================================================================
# Projected BGC zones - INSET
#===============================================================================

## parameters
grid <- "Salish1"
grid.dem <- "dem1_Salish"
rcp="rcp45"
proj.year=2055

# fplot=paste("InputData\\",grid, "_15GCM-Ensemble_rcp45_2055MSY.csv", sep="")
# 
# Y0 <- fread(fplot, select = Columns, stringsAsFactors = FALSE, data.table = FALSE) #fread is faster than read.csv
# 
# Y0 <- Y0[!is.na(Y0[,2]),]
# 
# #####generate some additional variables
# Y0$PPT_MJ <- Y0$PPT05 + Y0$PPT06 # MaY/June precip
# Y0$PPT_JAS <- Y0$PPT07 + Y0$PPT08 + Y0$PPT09 # July/Aug/Sept precip
# Y0$PPT.dormant <- Y0$PPT_at + Y0$PPT_wt # for calculating spring deficit
# Y0$CMD.def <- 500 - (Y0$PPT.dormant)# start of growing season deficit original value was 400 but 500 seems better
# Y0$CMD.def [Y0$CMD.def < 0] <- 0 #negative values set to zero = no deficit
# Y0$CMDMax <- Y0$CMD07
# Y0$CMD.total <- Y0$CMD.def + Y0$CMD
# 
# ##Predict future subzones######
#   BGC.pred.inset <- as.character(predict(BGCmodel, Y0))
# write.csv(BGC.pred.inset, paste("OutputData\\BGC.pred", grid, rcp, proj.year, "csv", sep="."), row.names = F)

BGC.pred.inset <- read.csv(paste("OutputData\\BGC.pred", grid, rcp, proj.year, "csv", sep="."))[,1]
unique(BGC.pred.inset)



######################
##reduce subzone-variant to zone

zone.pred.ref <- gsub("[:a-z:]","",BGC.pred.ref) 
zone.pred.ref <- gsub("[:1-9:]","",zone.pred.ref) 

zone.pred <- gsub("[:a-z:]","",BGC.pred) 
zone.pred <- gsub("[:1-9:]","",zone.pred) 

zone <- sort(unique(c(zone.pred.ref, zone.pred)))
zone <- factor(zone, levels=zone)

zone.pred.ref <- factor(zone.pred.ref, levels=levels(zone))
zone.pred <- factor(zone.pred, levels=levels(zone))

BGC.pred.ref.inset <- factor(as.character(BGC.pred.ref.inset), levels=levels(BGC.inset))
BGC.pred.inset <- factor(as.character(BGC.pred.inset), levels=levels(BGC.inset))

BGC.inset <- sort(unique(c(as.character(BGC.pred.ref.inset), as.character(BGC.pred.inset))))
BGC.inset <- factor(BGC.inset, levels=BGC.inset)







##############
# (A) reference BGC map
##############

#BGC zones
BGCcolors <- read.csv("C:\\Colin\\Projects\\2019_CCISS\\InputData\\BGCzone_Colorscheme.csv")
colors = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)][-1]

ColScheme.zone <- rep(NA, length(zone))
ColScheme.zone[which(zone%in%BGCcolors$zone)] <- as.character(BGCcolors$HEX[match(zone[which(zone%in%BGCcolors$zone)], BGCcolors$zone)])
set.seed(8)
ColScheme.zone[-which(zone%in%BGCcolors$zone)] <- sample(colors, length(zone[-which(zone%in%BGCcolors$zone)]))
ColScheme.zone <- factor(ColScheme.zone, levels=ColScheme.zone)

set.seed(1)
ColScheme.inset <- sample(colors, length(BGC.inset))
ColScheme.inset <- factor(ColScheme.inset, levels=ColScheme.inset)


grid <- "WNA2"
grid.dem <- "dem2_WNA"
grid.data <- read.csv(paste("InputData\\", grid, ".csv", sep = ""))
dem <- raster(paste("InputData\\", grid.dem,".tif", sep=""))
land.fine <- which(!is.na(values(dem)))  # raster cells with no dem value
P4S.latlon <- CRS("+proj=longlat +datum=WGS84")
X <- dem
values(X) <- NA

values(X) <- NA
values(X)[land.fine] <- zone.pred.ref
values(X)[1:length(levels(zone.pred.ref))] <- 1:length(levels(zone.pred.ref)) # this is a patch that is necessary to get the color scheme right. 

png(filename=paste("Results\\Manu_BGCprojections\\BGCzones", "png",sep="."), type="cairo", units="in", width=6.5, height=8, pointsize=9, res=600)

par(mar=c(0.1,0.1,0.1,0.1))
# image(hill, xlim=c(-135, -108), ylim=c(39, 60), col=alpha(grey(0:100/100), 1), xaxt="n", yaxt="n", maxpixels= ncell(hill))
image(X, xlim=c(-135, -108), ylim=c(39, 60), xaxt="n", yaxt="n", col=alpha(ColScheme.zone, 1), maxpixels=ncell(X))
plot(mask.NA, add=T, col="white", border=F)
plot(bdy.usa1, add=T, lwd=0.8)
plot(bdy.can1, add=T, lwd=0.8)

bgcs <- levels(zone.pred.ref)
for(bgc in bgcs){
  temp <- grid.data[which(zone.pred.ref==bgc),]
  median.lat <- temp$lat[max(which(temp$lat >= median(temp$lat)))]
  pt <- round(median(which(temp$lat==median.lat)))
  points(temp[pt,4:3], pch=21, bg=alpha(ColScheme.zone[which(bgcs==bgc)], 1), cex=1.5)
  text(temp[pt,4:3]-c(0,0), bgc, pos=4, cex=1, font=2)
  print(paste(which(bgcs==bgc), "-", bgc))
}
box() 
mtext(paste("(A) ", sep=""), side=3, line=-14.5, adj=0.02, cex=1.5, font=2)

grid <- "Salish1"
grid.dem <- "dem1_Salish"
grid.data <- read.csv(paste("InputData\\", grid, ".csv", sep = ""))
dem <- raster(paste("InputData\\", grid.dem,".tif", sep=""))
land.fine <- which(!is.na(values(dem)))  # raster cells with no dem value
length(land.fine)

extent <- c(-124.75, -121,46.5,50)
dem <- crop(dem, extent)
land.fine <- which(!is.na(values(dem)))  # raster cells with no dem value

select.crop <- which(grid.data$lon>extent[1] & grid.data$lon<extent[2] & grid.data$lat>extent[3] & grid.data$lat<extent[4])
grid.data.crop <- grid.data[select.crop, ]

X <- dem
values(X)[land.fine] <- BGC.pred.ref.inset[select.crop]
values(X)[1:length(levels(BGC.inset))] <- 1:length(levels(BGC.inset)) # this is a patch that is necessary to get the color scheme right. 
xlim=c(extent(X)[1], extent(X)[2])
ylim=c(extent(X)[3], extent(X)[4])
rect(xlim[1],  ylim[1],  xlim[2],  ylim[2],  col=(alpha("white", 0)), lwd=1.5)

par(plt = c(0.01, 0.375, 0.01, 0.45), new = TRUE)
image(X, xaxt="n", yaxt="n", xlim=xlim, ylim=ylim, col=alpha(ColScheme.inset, 1), maxpixels=ncell(X)) 
plot(mask.NA, add=T, col="white", border=F)
plot(bdy.usa1, add=T, lwd=0.8)
plot(bdy.can1, add=T, lwd=0.8)
mtext(paste("(B) ", sep=""), side=1, line=-1.5, adj=0.02, cex=1.5, font=2)

bgcs <- BGC.inset
for(bgc in bgcs){
  temp <- grid.data.crop[which(BGC.pred.ref.inset[select.crop]==bgc),]
  median.lat <- temp$lat[max(which(temp$lat >= median(temp$lat)))]
  pt <- round(quantile(which(temp$lat==median.lat), 0.2))
  points(temp[pt,4:3], pch=21, bg=alpha(ColScheme.inset[which(bgcs==bgc)], 1), cex=1.5)
  text(temp[pt,4:3]-c(0,0), bgc, pos=4, cex=.8, font=2)
  print(paste(which(bgcs==bgc), "-", bgc))
}

box()
dev.off()
