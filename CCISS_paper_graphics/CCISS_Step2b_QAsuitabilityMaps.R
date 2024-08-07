
##======================================================================================
## CCISS Publication Scripts
## Step 2b - Suitability maps for quality control of the suitability ratings
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


# Knowledge Tables
treesuit="TreeSpp_ESuit_v11_18"
SiteSeries_Use <-read.csv(paste("InputData/","SiteSeries_Use_5",".csv",sep=""),stringsAsFactors=FALSE,na.strings=".")
spps.lookup <- read.csv("InputData\\Tree speciesand codes_2.0_2May2019.csv")


# #==================================================
# # create the dem and ClimateWNA file from scratch
# #==================================================
# 
# setwd("C:\\Colin\\Projects\\2019_CCISS")
# 
# r <- raster("C:\\Colin\\SpatialData\\DEM\\namer_dem1.bil")
# plot(r, xlim=c(-140, -108), ylim=c(39,60))
# plot(bdy.usa, add=T)
# plot(bdy.can, add=T)
# r.crop <- crop(r, c(-140, -108,39,60))
# r.crop <- aggregate(r.crop, fact=2)
# writeRaster(r.crop, filename="InputData\\dem2_WNA.tif", format="GTiff", overwrite=TRUE)
# ## create a data frame of the projected coordinates of the dem
# r.pts <- rasterToPoints(r.crop, spatial=T) #create a spatial points data frame of the non-NA values of the DEM
# r.pts <- as.data.frame(r.pts) #projected coordinates of the dem
# 
# ## create the climateNA input file
# CNAinput <- data.frame(id1=1:dim(r.pts)[1], id2=rep(NA, dim(r.pts)[1]), lat=r.pts$y, lon=r.pts$x, el=r.pts$namer_dem1)
# str(CNAinput)
# write.csv(CNAinput,"InputData\\WNA2.csv", row.names=FALSE)


# setwd("C:\\Colin\\Projects\\2019_CCISS")
# r <- raster("C:\\Colin\\SpatialData\\DEM\\namer_dem1.bil")
# par(mar=c(2,2,0,0))
# plot(r, xlim=c(-124.75, -121), ylim=c(46,50.5))
# # plot(bdy.usa, add=T)
# # plot(bdy.can, add=T)
# r.crop <- crop(r, c(-124.75, -121,46,50.5))
# # r.crop <- aggregate(r.crop, fact=2)
# writeRaster(r.crop, filename="InputData\\dem1_Salish.tif", format="GTiff", overwrite=TRUE)
# ## create a data frame of the projected coordinates of the dem
# r.pts <- rasterToPoints(r.crop, spatial=T) #create a spatial points data frame of the non-NA values of the DEM
# r.pts <- as.data.frame(r.pts) #projected coordinates of the dem
# 
# ## create the climateNA input file
# CNAinput <- data.frame(id1=1:dim(r.pts)[1], id2=rep(NA, dim(r.pts)[1]), lat=r.pts$y, lon=r.pts$x, el=r.pts$namer_dem1)
# str(CNAinput)
# write.csv(CNAinput,"InputData\\Salish1.csv", row.names=FALSE)

#==================================================
# BGC projections (commented out because only have to do it once)
#==================================================
setwd("C:\\Colin\\Projects\\2019_CCISS")

## parameters
grid <- "WNA2"
grid.dem <- "dem2_WNA"

# ###Load random forest model
# fname="Rcode\\BGCv11_AB_USA_16VAR_SubZone_RFmodel.Rdata"
# # fname="Rcode\\BGCv11_AB_USA_500pts_27VAR_Subzone_RFmodel.Rdata"
# #fname = (file.choose())
# load(fname)
# 
# 
# fplot=paste("InputData\\", grid, "_Normal_1961_1990MSY.csv", sep="")
# 
# 
# Y0 <- fread(fplot, stringsAsFactors = FALSE, data.table = FALSE) #fread is faster than read.csv
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
# ##Predict reference subzones###### (need to break it up because it is big)
# BGC.pred.ref <- vector()
# factor <- 9
# for(i in 1:factor){
#   temp <- Y0[(dim(Y0)[1]/factor*(i-1)+1):(dim(Y0)[1]/factor*i),-c(1:3)]
#   # str(temp)
#   BGC.pred.ref <- c(BGC.pred.ref, as.character(predict(BGCmodel, temp)))
# print(i)
# }
# BGC.pred.ref
# write.csv(BGC.pred.ref, paste("OutputData\\BGC.pred.ref", grid, "csv", sep="."), row.names = F)

BGC.pred.ref <- read.csv(paste("OutputData\\BGC.pred.ref", grid, "csv", sep="."))[,1]
unique(BGC.pred.ref)

#==================================================
# spatial data
#==================================================
grid <- read.csv(paste("InputData\\", grid, ".csv", sep = ""))
dem <- raster(paste("InputData\\", grid.dem,".tif", sep=""))
land.fine <- which(!is.na(values(dem)))  # raster cells with no dem value
P4S.latlon <- CRS("+proj=longlat +datum=WGS84")
X <- dem
values(X) <- NA

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

### admin boundaries
bdy.bc <- readOGR("InputData\\BC_AB_US_Shp\\ProvincialOutline.shp")

#==================================================
# create a hillshade
#==================================================

r <- raster("C:\\Colin\\SpatialData\\DEM\\namer_dem1.bil")
r <- crop(r, c(-140, -108,39,60))
values(r)[is.na(values(r))] <- 0
slope = terrain(r, opt='slope')
aspect = terrain(r, opt='aspect')
hill = hillShade(slope, aspect, 40, 270)

#===============================================================================
# find the species suitability each projection/edatope/species combination
#===============================================================================

SiteLookup <- data.frame(MergedBGC=unique(SiteSeries_Use$MergedBGC))
edatopes<- c("B2", "C4", "D6")
for(edatope in edatopes){
  # SiteLookup <- cbind(SiteLookup, SiteSeries_Use$SS_NoSpace[match(SiteLookup[,1], SiteSeries_Use$MergedBGC[which(SiteSeries_Use$Use==edatope)])])
  SiteLookup <- cbind(SiteLookup, SiteSeries_Use$SS_NoSpace[which(SiteSeries_Use$Edatopic==edatope)])
  names(SiteLookup)[which(edatopes==edatope)+1] <- edatope
}
str(SiteLookup)

# Import suitability tables
S1 <- read.csv(paste("InputData/",treesuit,".csv",sep=""),stringsAsFactors=F,na.strings=".")
S1 <- unique(S1)

# ## EDA: are there suitabilities for all projected units? 
# NoSuit <- PredSum$Group.1[-which(PredSum$Group.1%in%S1$BGC)]
# NoSuit[rev(order(NoSuit$x)),]
# 
# ## EDA: Which site series are missing suitabilities? 
# for(edatope in edatopes) assign(paste("NoSuit", edatope, sep="."), SiteLookup[-which(SiteLookup[,which(names(SiteLookup)==edatope)]%in%S1$Unit),which(names(SiteLookup)==edatope)])
# for(edatope in edatopes) print(get(paste("NoSuit", edatope, sep=".")))

# select the species to run the analysis on
spps <- unique(S1$Spp)
spps <- spps[-which(spps=="X")]
spps.candidate <- spps.lookup$TreeCode[-which(spps.lookup$Exclude=="x")]
spps <- spps[which(spps%in%spps.candidate)] 

spp="Fd"
edatope="C4"
for(spp in spps){
  for(edatope in edatopes){
    # get the suitability for the reference period predicted BGC. 
    BGC.pred <- as.character(BGC.pred.ref) # get the BGC prediction 
    # BGC.pred[which(BGC.pred%in%Crosswalk$Modeled)] <- as.character(Crosswalk$Tables[match(BGC.pred[which(BGC.pred%in%Crosswalk$Modeled)], Crosswalk$Modeled)]) # XXX THIS IS NOT CORRECT. NEED TO FIGURE OUT HOW TO INCORPORATE THE CROSSWALK TABLE PROPERLY. sub in the crosswalk between the modeled units and the table units
    # get the suitability for the selected species associated with each site series
    suit <- S1$ESuit[which(S1$Spp==spp)][match(as.vector(unlist(SiteLookup[which(names(SiteLookup)==edatope)])), S1$Unit[which(S1$Spp==spp)])]
    # suit <- S1$Suitability[which(S1$Spp==spp)][match(as.vector(unlist(SiteLookup[which(names(SiteLookup)==edatope)])), S1$Unit[which(S1$Spp==spp)])]
    Suit.ref <- as.numeric(suit[match(BGC.pred, SiteLookup$MergedBGC)])
    Suit.ref[Suit.ref>3] <- NA
    assign(paste("Suit.ref", spp, edatope, "csv", sep="."), Suit.ref)
    # values(X)[land.fine] <- Suit.ref
    # writeRaster(X, paste("OutputData\\Suit.ref", grid, spp, edatope,"tif", sep="."),overwrite=TRUE, format="GTiff")
    
    print(edatope)
  }
  print(spp)
}


#===============================================================================
# map the suitability for each species
#===============================================================================
edatope.name <- c("Subxeric-poor", "Mesic-medium", "Hygric-rich")

breakseq <- c(0.5,1.5,2.5,3.5,5)
# ColScheme <- c(brewer.pal(9,"Greys")[9], brewer.pal(9,"Blues")[7], brewer.pal(9,"Greens")[4], "white")
ColScheme <- c("darkgreen", "dodgerblue1", "gold2", "white")
length(ColScheme)

for(spp in spps){
  for(edatope in edatopes){
    
    png(filename=paste("Results\\QAsuitability\\EnvSuit",spp, edatope, "png",sep="."), type="cairo", units="in", width=8.5, height=10, pointsize=10, res=300)
    
    values(X) <- NA
    values(X)[land.fine] <- get(paste("Suit.ref", spp, edatope, "csv", sep="."))

    par(mar=c(0.1,0.1,0.1,0.1))
    image(hill, xlim=c(-135, -108), ylim=c(39, 60), col=alpha(grey(0:100/100), 0.7), xaxt="n", yaxt="n", maxpixels= ncell(hill))
    image(X, add=T, xlim=c(-125, -112.5), ylim=c(40, 49), xaxt="n", yaxt="n", col=alpha(ColScheme, 0.6), breaks=breakseq , legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE)
    plot(mask.NA, add=T, col="white", border=F)
    plot(bdy.usa1, add=T, lwd=0.4)
    plot(bdy.can1, add=T, lwd=0.4)
    legend("bottomleft", legend=c("1 (primary)", "2 (secondary)", "3 (tertiary)"), fill=ColScheme, bty="n", cex=1.3, inset=0.1, title="Environmental\nsuitability")
    
    Common <- as.character(spps.lookup$EnglishName[which(spps.lookup$TreeCode==spp)])
    Latin <- as.character(spps.lookup$ScientificName[which(spps.lookup$TreeCode==spp)])
    text(-134.5, 46.5, spp, font=2, pos=4, cex=2)
    if(length(Common)>0){
    text(-134.5, 46, Common, font=2, pos=4, cex=1.6)
    text(-134.5, 45.5, Latin, font=3, pos=4, cex=1.3)
    }
    text(-134.5, 45, paste("Edatope: ", edatope, " (", edatope.name[which(edatopes==edatope)], ")", sep=""), font=1, pos=4, cex=1.3)

    box()
    
    dev.off()
    
    print(edatope)
  }
  print(spp)
}


