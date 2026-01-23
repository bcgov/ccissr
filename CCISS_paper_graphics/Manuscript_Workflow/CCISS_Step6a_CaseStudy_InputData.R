
## CCISS case study for a single BGC unit

## May 2019

## Colin Mahony 778-288-4008 c_mahony@alumni.ubc.ca


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

library(bcmaps)
data.frame(available_layers())

P4S.latlon <- CRS("+proj=longlat +datum=WGS84")
P4S.albers <- CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs")



# studyarea <- "Tahltan"
# studyarea <- "LakesTSA"
# studyarea <- "WetzinKwa"
# studyarea <- "Sunshine"
# studyarea <- "Kamloops"
# studyarea <- "Quesnel"
# studyarea <- "Cassiar"
# studyarea <- "FortStJohn"
# studyarea <- "PrinceGeorge"
# studyarea <- "JohnPrinceRF"
# studyarea <- "Cranbrook"
studyarea <- "Nimpkish"

# tsa <- tsa(class = "sp")
# unique(tsa$TSA_NUMBER_DESCRIPTION)


outdir <- paste("C:/Users/CMAHONY/OneDrive - Government of BC/Shiny_Apps/ccsummary-", studyarea, sep="")


# ## create a simplified BGC map for use in the locator
# library(bcmaps)
# # install.packages('bcmapsdata', repos='https://bcgov.github.io/drat/')
# library(bcmapsdata)
# library(rmapshaper)
# bec <- bec()
# bec.simpler <- rmapshaper::ms_simplify(input = as(bec, 'Spatial')) %>%
#   st_as_sf()
# bec.simpler <- st_transform(bec.simpler, CRS("+proj=longlat +datum=WGS84 +no_defs"))
# bec.simpler <- bec.simpler[,c(3,8)]
# st_write(bec.simpler, dsn="data/bec.simpler.shp", append=F)
# 
# ggplot() +
#   geom_sf(data = bec.simple[bec.simple$MAP_LABEL == "BGxh1",],
#           aes(fill = ZONE, col = ZONE)) +
#   scale_fill_manual(values = bec_colors()) +
#   scale_colour_manual(values = bec_colours())

#===============================================================================
# Set analysis Parameters
#===============================================================================

source("./_CCISS_Packages.R") ## packages required
source("./_CCISS_Functions.R") ## common functions
source("./_CCISS_Parameters.R") ## settings used through all scripts


rcp.focal="rcp45"
proj.year.focal=2055

#BGC zone color scheme
BGCcolors$colour <- as.character(BGCcolors$colour)
BGCcolors$colour[match(BGCcolors.BC$zone, BGCcolors$classification)] <- as.character(BGCcolors.BC$HEX)
ColScheme <- factor(BGCcolors$colour, levels=BGCcolors$colour)
zones <- factor(BGCcolors$classification, levels=BGCcolors$classification)

GCMs <-  c("CanESM2","CCSM4","CESM1-CAM5","CSIRO-Mk3-6-0", "GISS-E2R", "HadGEM2-ES", "MIROC-ESM", "MPI-ESM-LR")    

### Load random forest model
load(fname)
#rownames(importance(BGCmodel2)) ### shows the variable used in the RFmodel
vars <- as.data.frame(BGCmodel2$variable.importance)
vars <- row.names(vars)
# setwd('C:/GitHub/2019_CCISS')

Columns <- unique(c("PPT05", "PPT06", "PPT07", "PPT08", "PPT09", "PPT_at", 
                    "PPT_wt", "CMD07", "CMD", "MAT", "PPT_sm", "Tmin_wt", "Tmax_sm",
                    vars[!vars %in% c("PPT_MJ", "PPT_JAS", "PPT.dormant", "CMD.def", "CMDMax", "CMD.total")]))


#===============================================================================
# DEM and boundary
#===============================================================================

### lakes
lakes <- readOGR("C:/Users/CMAHONY/OneDrive - Government of BC\\SpatialData\\50k_layers\\LakesBig.shp")
lakes <- gSimplify(lakes, 10)
# projection(lakes) <- P4S.AEA

temp <- as.data.frame(available_layers())
## boundary for the study area
# bdy.aea <- readOGR("C:/GitHub/2019_CCISS/inputs/shapes/Tahltan/TahltanLUP_BoundaryV7_200702.shp")
# bdy.aea <- readOGR("C:/GitHub/2019_CCISS/inputs/shapes/Lakes/LakesTSA_3L3P_190815.shp")
# bdy.aea <- readOGR("C:/GitHub/2019_CCISS/inputs/shapes/WetzinKwa/UBM_LUP.shp")
# bdy.aea <- readOGR("C:/Users/CMAHONY/OneDrive - Government of BC/SpatialData/Boundaries/WL_ScheduleB_all.shp")
# bdy.aea <- tsa[grep(substr(studyarea, 1,6), tsa$TSA_NUMBER_DESCRIPTION),]
bdy.aea <- readOGR(paste(outdir, "bdy/Nimpkish_LSU.shp", sep="/"))
# bdy.aea <- readOGR(paste(outdir, "bdy/bdy.shp", sep="/"))

bdy.aea <- gSimplify(bdy.aea, 10)
bdy.aea <- gUnaryUnion(bdy.aea)
bc_bound <- bc_bound_hres(class="sp")
plot(bdy.aea)
plot(bc_bound, add=T, border=2)
plot(lakes, add=T, border=2)
bdy.aea <- gIntersection(bdy.aea,bc_bound)

bdy <- spTransform(bdy.aea, P4S.latlon) #reproject the gridpoints
plot(bdy)

# # ####### create a polygon mask for the focal unit 
# border=0.5
# my_box = as(extent(buffer(bdy, border)), "SpatialPolygons")      		# convert extent box to shapefile (rectangle)
# proj4string(my_box) = projection(bdy)				# assign spatial projection to extent object
# focal.mask <- gDifference(my_box, bdy)
# projection(focal.mask)  #verify latlong projection of the study area boundary
# par(mar=c(0,0,0,0))
# plot(focal.mask)
# plot(my_box)
# plot(bdy, add=T)

plot(bdy)
# bdy.simplify <- gSimplify(bdy, 0.001)
# bdy.simplify <- gUnaryUnion(bdy.simplify)
plot(bdy)
writeOGR(as(bdy, "SpatialPolygonsDataFrame" ), dsn = paste(outdir, "/bdy/bdy", studyarea, "shp", sep="."), layer="TSA_NUMBER", driver = "ESRI Shapefile")

## create the study area

library(rcanvec)
nts.bybbox(sp::bbox(bdy), "nts.SCALESERIES")
mapsheets <- ntsstring(bbox=sp::bbox(bdy))
mapsheets <- unique(substr(mapsheets,1,4))
mapsheets <- paste(mapsheets, collapse="|")
dir <- "C:\\Users/CMAHONY/OneDrive - Government of BC\\SpatialData\\CDEM/unzipped"
files <- list.files(dir, pattern="*.tif")
files <- files[grep(mapsheets, files)]
for(i in 1:length(files)){
  temp <- aggregate(raster(paste(dir, files[i], sep = "/")), fact=8) #reduce resolution
  if(i==1) alt <- temp else alt <- mosaic(alt, temp, fun=mean)
  print(i)
}
alt <- crop(alt,bdy)
par(mar=c(0,0,0,0))
plot(alt)
plot(bdy, add=T)
r <- mask(alt, bdy)
plot(r)
alt

par(mar=c(0,0,0,0))
    plot(r)
plot(bdy, add=T)

# r <- raster(paste("inputs\\dem1_", studyarea, ".tif", sep=""))
writeRaster(r, filename=paste("inputs\\dem1_", studyarea, ".tif", sep=""), format="GTiff", overwrite=TRUE)
## create a data frame of the projected coordinates of the dem
r.pts <- rasterToPoints(r, spatial=T) #create a spatial points data frame of the non-NA values of the DEM
projection(r.pts) <- P4S.latlon 
str(r.pts)

### Attribute the points with BGC linework
BGC <- bec(class="sp")
BGC <- crop(BGC, buffer(bdy.aea, 4000))

pts.aea <- spTransform(r.pts, P4S.albers) # reproject to lat-long 
# plot(pts.aea, add=T, pch=16, cex=0.5)
BGC.pts <- over(pts.aea,BGC)  ##extract the BGCv11 attributes to the points (this takes a while, like 30ish minutes)
str(BGC.pts)

r.pts <- as.data.frame(r.pts) #projected coordinates of the dem
instudyarea <- unlist(cellFromPolygon(r, bdy))
nonland <- which(is.na(values(r)[instudyarea]))
land <- if(length(nonland)>0) instudyarea[-which(is.na(values(r)[instudyarea]))] else instudyarea
length(instudyarea)
length(land)

## create the climateNA input file
CNAinput <- data.frame(id1=land, id2=as.character(BGC.pts$MAP_LABEL), lat=r.pts$y, lon=r.pts$x, el=r.pts[,1])
str(CNAinput)
write.csv(CNAinput,paste("inputs\\", studyarea, ".csv", sep=""), row.names=FALSE)
#NOTE! FUTURE NORMALS SHOULD BE QUERIED FROM CLIMATEBC AS BATCH FILES FOR ALL GCMS AND FOR ONLY ONE RCP AND NORMAL PERIOD. 



# #===============================================================================
# # create a hillshade backdrop
# #===============================================================================
# 
# 
# dem <- raster(paste("inputs\\dem1_",studyarea,".tif", sep=""))
# 
# library(rcanvec)
# nts.bybbox(sp::bbox(bdy), "nts.SCALESERIES")
# mapsheets <- ntsstring(bbox=sp::bbox(bdy))
# mapsheets <- unique(substr(mapsheets,1,4))
# mapsheets <- paste(mapsheets, collapse="|")
# dir <- "C:\\Users/CMAHONY/OneDrive - Government of BC\\SpatialData\\CDEM/unzipped"
# files <- list.files(dir, pattern="*.tif")
# files <- files[grep(mapsheets, files)]
# for(i in 1:length(files)){
#   temp <- aggregate(raster(paste(dir, files[i], sep = "/")), fact=8) #reduce resolution
#   if(i==1) alt <- temp else alt <- mosaic(alt, temp, fun=mean)
#   print(i)
# }
# # alt <- crop(alt,X)
# slope = terrain(alt, opt='slope')
# aspect = terrain(alt, opt='aspect')
# hill = hillShade(slope, aspect, 40, 270)
# plot(hill, col=grey(0:100/100), legend=FALSE)
# 
# ## subset the raster for focal cells only
# focal.buffer <- gBuffer(focal.mask, width=-res(dem)[1])
# focal.cells.buffer <- which(is.na(values(mask(dem,focal.buffer))))
# focal.cells <- which(is.na(values(mask(dem,focal.mask))))
# 
# X <- dem
# values(X) <- NA
# values(X)[focal.cells] <- 1
# plot(X)
# 
# plot(dem)
# plot(bdy, add=T, border="blue")
# plot(hill, col=grey(0:100/100), legend=FALSE)
# plot(focal.mask, add=T, col=alpha("white", 0.5))

#===============================================================================
#===============================================================================
# BGC Projections
#===============================================================================
#===============================================================================


bdy <- readOGR(paste(outdir, "/bdy/bdy", studyarea, "shp", sep="."))
grid <- read.csv(paste("inputs\\", studyarea, ".csv", sep=""), stringsAsFactors=F)
dem <- raster(paste("inputs\\dem1_",studyarea,".tif", sep=""))
# land <- grid$id1
land <- which(!is.na(values(dem)))
# instudyarea <- unlist(cellFromPolygon(dem, bdy))
# nonland <- which(is.na(values(dem)[instudyarea]))
# land <- if(length(nonland)>0) instudyarea[-which(is.na(values(dem)[instudyarea]))] else instudyarea
X <- dem
plot(X)
values(X) <- NA

bgc <- grid[,2]
zone <- rep(NA, length(bgc))
for(i in BGCcolors.BC$zone){ zone[grep(i,bgc)] <- i }
table(zone)

length(bgc)
length(land)

values(X)[land] <- factor(bgc)
plot(X)

#===============================================================================
# create a reference BGC map
#===============================================================================


png(filename=paste(outdir, "/www/refmap", studyarea,"zones.png",sep="."), type="cairo", units="in", width=4, height=4, pointsize=11, res=200)
par(mar=c(0,0,0,0))
values(X)[land] <-  factor(zone, levels=zones)
values(X)[1:length(zones)] <- 1:length(zones) # this is a patch that is necessary to get the color scheme right.
plot(X, xaxt="n", yaxt="n", col=alpha(ColScheme, 1), legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE) 
values(X)[-(1:length(zones))] <- NA # cover up the color bar
image(X, add=T, col="white") # cover up the color bar
plot(bdy, add=T, lwd=1)
# box()
dev.off()

png(filename=paste(outdir, "/www/refmap", studyarea,"variants.png",sep="."), type="cairo", units="in", width=4, height=4, pointsize=11, res=200)
par(mar=c(0,0,0,0))
values(X)[land] <-  factor(bgc, levels=BGCcolors.subzone$classification)
values(X)[1:length(BGCcolors.subzone$classification)] <- 1:length(BGCcolors.subzone$classification) # this is a patch that is necessary to get the color scheme right.
plot(X, xaxt="n", yaxt="n", col=alpha(BGCcolors.subzone$colour, 1), legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE) 
values(X)[-(1:length(BGCcolors.subzone$classification))] <- NA # cover up the color bar
image(X, add=T, col="white") # cover up the color bar
plot(bdy, add=T, lwd=1) 
# box()
dev.off()


#===============================================================================
# BGC Projections for reference period
#===============================================================================

# setwd("C:\\Users/CMAHONY/OneDrive - Government of BC\\Projects\\2019_CCISS")

fplot=paste("inputs\\", studyarea, "_Normal_1961_1990MSY.csv", sep="")

Y0 <- fread(fplot, stringsAsFactors = FALSE, data.table = FALSE) #fread is faster than read.csv
#####generate some additional variables
Y0$PPT_MJ <- Y0$PPT05 + Y0$PPT06 # MaY/June precip
Y0$PPT_JAS <- Y0$PPT07 + Y0$PPT08 + Y0$PPT09 # July/Aug/Sept precip
Y0$PPT.dormant <- Y0$PPT_at + Y0$PPT_wt # for calculating spring deficit
Y0$CMD.def <- 500 - (Y0$PPT.dormant)# start of growing season deficit original value was 400 but 500 seems better
Y0$CMD.def [Y0$CMD.def < 0] <- 0 #negative values set to zero = no deficit
Y0$CMDMax <- Y0$CMD07
Y0$CMD.total <- Y0$CMD.def + Y0$CMD

##Predict future subzones######
BGC.pred.ref <- predict(BGCmodel2, Y0)
BGC.pred.ref <- BGC.pred.ref$predictions

values(X)[land] <- BGC.pred.ref
writeRaster(X, paste(outdir, "/data/BGC.pred", studyarea, "ref.tif", sep="."),overwrite=TRUE, format="GTiff")

## calculate mean climate of study area for use in calculating change
Y.refmean <- apply(Y0[,-c(1:5)], 2, FUN=mean, na.rm=T)
write.csv(t(as.data.frame(Y.refmean)), paste(outdir, "/data/clim.refMean",studyarea,"csv", sep="."), row.names = F)

change <- data.frame("GCM"="obs", "rcp"="obs", "proj.year"=1975, as.data.frame(t(rep(0, length(Y.refmean)))))
names(change)[-c(1:3)] <- names(Y.refmean)

#===============================================================================
# BGC Projections for recent period
#===============================================================================

# setwd("C:\\Users/CMAHONY/OneDrive - Government of BC\\Projects\\2019_CCISS")
hist.years <- c(1995, 2005, 2015)
hist.periods <- c("1991_2000", "2001_2010", "2011_2019")

# read in the data for the decades
for(hist.year in hist.years){
  hist.period <- hist.periods[which(hist.years==hist.year)]
  fplot=paste("inputs\\", studyarea, "_Decade_", hist.period, "MSY.csv", sep="")

  Y0 <- fread(fplot, stringsAsFactors = FALSE, data.table = FALSE)[-c(1:5)] #fread is faster than read.csv

  #####generate some additional variables
  Y0$PPT_MJ <- Y0$PPT05 + Y0$PPT06 # MaY/June precip
  Y0$PPT_JAS <- Y0$PPT07 + Y0$PPT08 + Y0$PPT09 # July/Aug/Sept precip
  Y0$PPT.dormant <- Y0$PPT_at + Y0$PPT_wt # for calculating spring deficit
  Y0$CMD.def <- 500 - (Y0$PPT.dormant)# start of growing season deficit original value was 400 but 500 seems better
  Y0$CMD.def [Y0$CMD.def < 0] <- 0 #negative values set to zero = no deficit
  Y0$CMDMax <- Y0$CMD07
  Y0$CMD.total <- Y0$CMD.def + Y0$CMD
  assign(paste("Y",hist.year,sep="."), Y0)
  print(hist.year)
}
# str(Y.2005)

## take the weighted mean of the decades
Y.hist <- (Y.1995*10+Y.2005*10+Y.2015*9)/29

##Predict future subzones######
BGC.pred.hist <- predict(BGCmodel2, Y.hist)
# write.csv(list(BGC.pred.hist$predictions), paste(outdir, "/data/BGC.pred", studyarea, "2005.csv", sep="."), row.names = F)
BGC.pred.hist <- BGC.pred.hist$predictions
values(X)[land] <- BGC.pred.hist
writeRaster(X, paste(outdir, "/data/BGC.pred", studyarea, "2005.tif", sep="."),overwrite=TRUE, format="GTiff")

## calculate climate change
Y.mean <- apply(Y.hist, 2, FUN=mean, na.rm=T)
change.temp <- Y.mean - Y.refmean
change <- rbind(change, data.frame("GCM"="obs", "rcp"="obs", "proj.year"=2005, as.data.frame(t(change.temp))))



#===============================================================================
# BGC Projections for future periods
#===============================================================================
rcp="rcp45"

# for(rcp in rcps){
  for(proj.year in proj.years){
    fplot=paste("inputs\\", studyarea, "_", rcp, "_", proj.year, "_MSY.csv", sep="")
    Y0 <- fread(fplot, stringsAsFactors = FALSE, data.table = FALSE) #fread is faster than read.csv
    
    #####generate some additional variables
    Y0$PPT_MJ <- Y0$PPT05 + Y0$PPT06 # MaY/June precip
    Y0$PPT_JAS <- Y0$PPT07 + Y0$PPT08 + Y0$PPT09 # July/Aug/Sept precip
    Y0$PPT.dormant <- Y0$PPT_at + Y0$PPT_wt # for calculating spring deficit
    Y0$CMD.def <- 500 - (Y0$PPT.dormant)# start of growing season deficit original value was 400 but 500 seems better
    Y0$CMD.def [Y0$CMD.def < 0] <- 0 #negative values set to zero = no deficit
    Y0$CMDMax <- Y0$CMD07
    Y0$CMD.total <- Y0$CMD.def + Y0$CMD
    
    ## assign single vectors to RCPs and proj.years
    Y0 <- separate(Y0, Year, into = c("GCM","rcp","proj.year"), sep = "_", remove = T)
    Y0$proj.year <- gsub(".gcm","",Y0$proj.year)
    
    GCMs <- unique(Y0[,1])
    for(GCM in GCMs){
      temp <- Y0[which(Y0[,1]==GCM & Y0[,2]==rcp & Y0[,3]==proj.year),]
      BGC.pred.proj <- predict(BGCmodel2, temp)
      # write.csv(list(BGC.pred.proj$predictions), paste(outdir, "/data/BGC.pred",studyarea, GCM, rcp, proj.year,"csv", sep="."), row.names = F)
      BGC.pred.proj <- BGC.pred.proj$predictions
      assign(paste("BGC.pred",GCM, rcp, proj.year,sep="."), BGC.pred.proj)
      
      ## calculate climate change
      Y.mean <- apply(temp[,which(names(temp)%in%names(Y.refmean))], 2, FUN=mean, na.rm=T)

      change.temp <- Y.mean - Y.refmean
      change <- rbind(change, data.frame("GCM"=GCM, "rcp"=rcp, "proj.year"=proj.year, as.data.frame(t(change.temp))))

      values(X)[land] <- BGC.pred.proj
      writeRaster(X, paste(outdir, "/data/BGC.pred", studyarea, GCM, rcp, proj.year,"tif", sep="."),overwrite=TRUE, format="GTiff")
      
      print(GCM)
    }
    print(proj.year)
  }
  # print(rcp)
# }

#calculate ensemble mean change
change.ensemble <- aggregate(change[-c(1,2),-c(1,2,3)], by=list(change[-c(1,2),2],change[-c(1,2),3]), FUN=mean, na.rm=T)
names(change.ensemble)[1:2] <- names(change)[2:3]
change <- rbind(change, data.frame("GCM"="ensemble", change.ensemble))

write.csv(change, paste(outdir, "/data/clim.meanChange",studyarea,"csv", sep="."), row.names = F)

#===============================================================================
# Ensemble Mean BGC Projections for future periods
#===============================================================================

BGC.pred <- BGC.pred.proj
levels.bgc <- levels(BGC.pred.proj)
write.csv(levels.bgc, paste(outdir, "/data/levels.bgc.csv", sep="."), row.names = F)

# determine vote winner BGC and ensemble agreement (WARNING: takes about 1 minute per rcp/proj.year)
# for(rcp in rcps){
  # proj.year=proj.years[2]
  for(proj.year in proj.years){
    temp <- as.data.frame(matrix(rep(NA, length(BGC.pred)*length(GCMs)), nrow=length(BGC.pred), ncol=length(GCMs)))
    for(GCM in GCMs){
      BGC.pred <- get(paste("BGC.pred",GCM, rcp, proj.year,sep="."))
      #add votes to votes matrix
      temp[,which(GCMs==GCM)] <- BGC.pred
      # print(GCM)
    }
    vote.winner <- function(x){return(names(which(table(x)==max(table(x))))[1])}
    agreement <- function(x){return(max(table(x)))}
    assign(paste("BGC.pred.ensemble", rcp, proj.year, sep="."), apply(temp, 1, vote.winner))
    # assign(paste("BGC.pred.agreement", rcp, proj.year, sep="."), apply(temp, 1, agreement))
    # write.csv(get(paste("BGC.pred.ensemble", rcp, proj.year, sep=".")), paste(outdir, "/data/BGC.pred",studyarea, "ensemble", rcp, proj.year,"csv", sep="."), row.names = F)
    # write.csv(agreement, paste(outdir, "/data/BGC.pred",studyarea, "agreement", rcp, proj.year,"csv", sep="."), row.names = F)

    values(X)[land] <- factor(get(paste("BGC.pred.ensemble", rcp, proj.year, sep=".")), levels=levels.bgc)
    writeRaster(X, paste(outdir, "/data/BGC.pred", studyarea, "ensemble", rcp, proj.year,"tif", sep="."),overwrite=TRUE, format="GTiff")

    print(proj.year)
  }
  # print(rcp)
# }


#===============================================================================
# Make and export summary tables of BGC units for each future
#===============================================================================

PredSum <- data.frame("GCM"="obs", "rcp"="obs", "proj.year"=1975, as.data.frame(table(BGC.pred.ref)))
names(PredSum)[4] <- "BGC.pred"

# Historical BGC
hist.years=2005
hist.year=2005
for(hist.year in hist.years){
  BGC.pred <- BGC.pred.hist
  assign(paste("BGC.pred", hist.year, sep="."), BGC.pred) #bgc projection
  PredSum <- rbind(PredSum, data.frame("GCM"="obs", "rcp"="obs", "proj.year"=2005, as.data.frame(table(BGC.pred))))
  print(hist.year)
}

# Future BGC
# for(rcp in rcps){
for(proj.year in proj.years){
  for(GCM in c("ensemble", GCMs)){
    BGC.pred <- get(paste("BGC.pred",GCM, rcp, proj.year,sep="."))
    PredSum <- rbind(PredSum, data.frame("GCM"=GCM, "rcp"=rcp, "proj.year"=proj.year, as.data.frame(table(BGC.pred))))
    # print(GCM)
  }
  print(proj.year)
}
# print(rcp)
# }

# write out summary of BGC units for each future. 
PredSum.wide <- reshape(PredSum, idvar = c(names(PredSum)[1:3]), timevar = "BGC.pred", direction = "wide")
write.csv(PredSum.wide, paste(outdir, "/data/PredSum.BGC",studyarea,"csv", sep="."), row.names = F)

#===============================================================================
#===============================================================================
# Species Feasibility Projections
#===============================================================================
#===============================================================================

# #===============================================================================
# # generate the vector of mapped BGCs
# #===============================================================================
 
points <- fread(paste("./inputs/", studyarea,".csv", sep=""))
BGC <- points$id2
BGC <- gsub(" ","",BGC)
sort(table(BGC))

#===============================================================================
# import will's lookup table for the site series associated with the selected edatope in each BGC unit
#===============================================================================

SiteLookup <- data.frame(BGC=unique(SiteSeries_Use$BGC))
edatopes<- c("B2", "C4", "D6")
for(edatope in edatopes){
  # SiteLookup <- cbind(SiteLookup, SiteSeries_Use$SS_NoSpace[match(SiteLookup[,1], SiteSeries_Use$MergedBGC[which(SiteSeries_Use$Use==edatope)])])
  SiteLookup <- cbind(SiteLookup, SiteSeries_Use$SS_NoSpace[which(SiteSeries_Use$Use==edatope)])
  names(SiteLookup)[which(edatopes==edatope)+1] <- edatope
}
str(SiteLookup)
write.csv(SiteLookup, "outputs/SiteLookup.csv", row.names = F)

#===============================================================================
# find the species suitability each projection/edatope/species combination
#===============================================================================

hist.years=2005

# Import suitability tables
S1 <- treesuit
S1 <- unique(S1)[,1:5]
dim(S1)
S1 <- unique(S1)
dim(S1)
write.csv(S1, "outputs/SuitLookup.csv", row.names = F)


## EDA: are there suitabilities for all projected units? 
NoSuit <- PredSum[-which(PredSum$BGC%in%S1$BGC),]
NoSuit[rev(order(NoSuit$Freq)),]

## EDA: Which site series are missing suitabilities? 
for(edatope in edatopes) assign(paste("NoSuit", edatope, sep="."), SiteLookup[-which(SiteLookup[,which(names(SiteLookup)==edatope)]%in%S1$SS_NoSpace),which(names(SiteLookup)==edatope)])
for(edatope in edatopes) print(get(paste("NoSuit", edatope, sep=".")))

## EDA: are there any units missing from the SiteSeries_Use table? 
BGClist <- unique(S1$BGC)
BGClist[-which(BGClist%in%SiteLookup$BGC)]

# select the species to run the analysis on
spps <- unique(S1$Spp)
spps <- spps[-which(spps=="X")]
spps.candidate <- spps.lookup$TreeCode[-which(spps.lookup$Exclude=="x")]
spps <- as.character(spps[which(spps%in%spps.candidate)] )

# for(spp in c("Pl", "Fd", "Cw", "Sx")){
for(spp in spps){
  for(edatope in edatopes){
    # get the suitability for the reference period predicted BGC.
    BGC.pred <- BGC.pred.ref
    
    # get the suitability for the selected species associated with each site series
    suit <- S1$Feasible[which(S1$Spp==spp)][match(as.vector(unlist(SiteLookup[which(names(SiteLookup)==edatope)])), as.character(S1$SS_NoSpace[which(S1$Spp==spp)]))]
    Suit.ref <- suit[match(BGC.pred, SiteLookup$BGC)]
    Suit.ref[is.na(Suit.ref)] <- 5 #set the NA values to suitability 5
    Suit.ref[Suit.ref==4] <- 5 #set 4 to suitability 5
    assign(paste("Suit.ref", spp, edatope, sep="."), Suit.ref)
    # write.csv(Suit.ref, paste(outdir, "/data/Suit.ref", studyarea, spp, edatope, "csv", sep="."), row.names = F)
    
    for(hist.year in hist.years){
      BGC.pred <- get(paste("BGC.pred", hist.year, sep="."))
      ## identify cells with no suitability interpretations
      bgc.exotic <- (1:length(BGC.pred))[-which(BGC.pred%in%unique(BGC))]
      bgc.exotic.noSuit <- bgc.exotic[-which(BGC.pred[bgc.exotic]%in%unique(S1$BGC))]
      
      suit <- S1$Feasible[which(S1$Spp==spp)][match(as.vector(unlist(SiteLookup[which(names(SiteLookup)==edatope)])), S1$SS_NoSpace[which(S1$Spp==spp)])]
      temp <- suit[match(BGC.pred, SiteLookup$BGC)]
      temp[is.na(temp)] <- 5 #set the NA values to suitability 5
      temp[temp==4] <- 5 #set 4 to suitability 5
      temp[bgc.exotic.noSuit] <- NA # set cells with no suitabilty interpretatoin to NA
      assign(paste("Suit", hist.year, spp, edatope, sep="."), temp)
      # write.csv(temp, paste(outdir, "/data/Suit", studyarea, hist.year, spp, edatope, "csv", sep="."), row.names = F)
      # print(hist.year)
    }
    
    # get the suitability for future periods, for each projection/edatope/species combination
    for(GCM in c(GCMs, "ensemble")){
      # for(rcp in rcps){
        for(proj.year in proj.years){
          # get the BGC projection and sub in the crosswalk between the modeled units and the table units
          BGC.pred <- get(paste("BGC.pred",GCM, rcp, proj.year,sep="."))
          # BGC.pred[which(BGC.pred%in%Crosswalk$Modeled)] <- as.character(Crosswalk$Tables[match(BGC.pred[which(BGC.pred%in%Crosswalk$Modeled)], Crosswalk$Modeled)]) # XXX THIS IS NOT CORRECT. NEED TO FIGURE OUT HOW TO INCORPORATE THE CROSSWALK TABLE PROPERLY. sub in the crosswalk between the modeled units and the table units

          ## identify cells with no suitability interpretations
          bgc.exotic <- (1:length(BGC.pred))[-which(BGC.pred%in%unique(BGC))]
          bgc.exotic.noSuit <- bgc.exotic[-which(BGC.pred[bgc.exotic]%in%unique(S1$BGC))]

          # get the suitability for the selected species associated with each site series
          suit <- S1$Feasible[which(S1$Spp==spp)][match(as.vector(unlist(SiteLookup[which(names(SiteLookup)==edatope)])), S1$SS_NoSpace[which(S1$Spp==spp)])]
          temp <- suit[match(BGC.pred, SiteLookup$BGC)]
          temp[is.na(temp)] <- 5 #set the NA values to suitability 5 (weights unsuitable a bit more heavily than suitable classes during averaging)
          temp[temp==4] <- 5 #set 4 to suitability 5
          temp[bgc.exotic.noSuit] <- NA # set cells with no suitabilty interpretation to NA
          assign(paste("Suit", GCM, rcp, proj.year, spp, edatope, sep="."), temp)
          # write.csv(temp, paste(outdir, "/delete/Suit", studyarea, GCM, rcp, proj.year, spp, edatope, "csv", sep="."), row.names = F)
          # print(proj.year)
        }
        # print(rcp)
      # }
      # print(GCM)
    }
    # print(edatope)
  }
  print(paste(spp, " (", round(which(spps==spp)/length(spps)*100, 0), "%)", sep=""))
}


#===============================================================================
# summarize the suitability of species for each scenario. 
#===============================================================================

## non-THLB BGCs for exclusion from results
points <- read.csv(paste("inputs\\",studyarea,".csv", sep=""))
BGC <- points$id2
BGC <- gsub(" ","",BGC)
exclude <- which(BGC%in%BGCs_notin_THLB$BGC[which(BGCs_notin_THLB$Exlude=="x")])
include <- if(length(exclude>0)) seq(1,length(BGC))[-exclude] else seq(1,length(BGC))

for(edatope in edatopes){
  
  #initiate tables to store summary values
  PredSum.suit <- data.frame(PredSum.wide[1:3], as.data.frame(matrix(rep(NA, length(spps)*dim(PredSum.wide)[1]), dim(PredSum.wide)[1])))
  names(PredSum.suit)[-c(1:3)] <- spps
  PredSum.spp <- PredSum.suit
  PredSum.suit.home <- PredSum.suit #home is for counting cells within historical range. 
  PredSum.spp.home <- PredSum.suit
  
  for(spp in spps){
    
    #reference period suitabilities
    Suit.ref <- get(paste("Suit.ref", spp, edatope, sep="."))[include]
    # Suit.ref <- read.csv(paste(outdir, "/data/Suit.ref", spp, edatope, "csv", sep="."))[include,1]
    # Suit.ref <- read.csv(paste(outdir, "/data/Suit.ref", studyarea, spp, edatope, "csv", sep="."))[,1]
    Suit.ref[Suit.ref==5] <- NA
    outRange.ref <- is.na(Suit.ref)
    Suit.ref[is.na(Suit.ref)] <- 5
    Suit.ref[Suit.ref==4] <- 3 #added this based on email from Will May 18, 2021
    Suit.ref <- 1-(Suit.ref-1)/4

    row <- 1
    col <- which(names(PredSum.suit)==spp)
    PredSum.suit[row,col] <- round(sum(Suit.ref))
    PredSum.spp[row,col] <- round(sum(Suit.ref>0))
    PredSum.suit.home[row,col] <- round(sum(Suit.ref[outRange.ref==F]))
    PredSum.spp.home[row,col] <- round(sum((Suit.ref>0)[outRange.ref==F]))
    
    for(hist.year in hist.years){
      Suit.proj <- get(paste("Suit", hist.year, spp, edatope, sep="."))[include]
      Suit.proj[is.na(Suit.proj)] <- 5
      Suit.proj[Suit.proj==4] <- 3 #added this based on email from Will May 18, 2021
      Suit.proj <- 1-(Suit.proj-1)/4

      row <- which(PredSum.suit$proj.year==hist.year)
      PredSum.suit[row,col] <- round(sum(Suit.proj))
      PredSum.spp[row,col] <- round(sum(Suit.proj>0))
      PredSum.suit.home[row,col] <- round(sum(Suit.proj[outRange.ref==F]))
      PredSum.spp.home[row,col] <- round(sum((Suit.proj>0)[outRange.ref==F]))
    }

    for(GCM in c(GCMs, "ensemble")){
      # for(rcp in rcps){
        for(proj.year in proj.years){
          Suit.proj <- get(paste("Suit", GCM, rcp, proj.year, spp, edatope, sep="."))[include]
          # Suit.proj <- read.csv(paste(outdir, "/delete/Suit", studyarea, GCM, rcp, proj.year, spp, edatope, "csv", sep="."))[,1][include]
          Suit.proj[is.na(Suit.proj)] <- 5
          Suit.proj[Suit.proj==4] <- 3 #added this based on email from Will May 18, 2021
          Suit.proj <- 1-(Suit.proj-1)/4
          
          row <- which(PredSum.suit$GCM==GCM & PredSum.suit$rcp==rcp & PredSum.suit$proj.year==proj.year)
          PredSum.suit[row,col] <- round(sum(Suit.proj))
          PredSum.spp[row,col] <- round(sum(Suit.proj>0))
          PredSum.suit.home[row,col] <- round(sum(Suit.proj[outRange.ref==F]))
          PredSum.spp.home[row,col] <- round(sum((Suit.proj>0)[outRange.ref==F]))
          # print(proj.year)
        }
        # print(rcp)
      # }
      # print(GCM)
    }
  print(paste(spp, " (", round(which(spps==spp)/length(spps)*100, 0), "%)", sep=""))
}

  assign(paste("PredSum.suit", edatope, sep="."), PredSum.suit)
  
  write.csv(PredSum.suit, paste(outdir, "/data/PredSum.suit",studyarea, edatope,"csv", sep="."), row.names = F)
  write.csv(PredSum.spp, paste(outdir, "/data/PredSum.spp",studyarea, edatope,"csv", sep="."), row.names = F)
  write.csv(PredSum.suit.home, paste(outdir, "/data/PredSum.suit.home",studyarea, edatope,"csv", sep="."), row.names = F)
  write.csv(PredSum.spp.home, paste(outdir, "/data/PredSum.spp.home",studyarea, edatope,"csv", sep="."), row.names = F)
  
    print(edatope)
}

#===============================================================================
# Write rasters of mean feasibilty change and binary appearance
#===============================================================================
edatope="C4"
for(edatope in edatopes){
  
# exclude insignificant species
suit.area <-  read.csv(paste(outdir, "/data/PredSum.suit",studyarea, edatope,"csv", sep="."))[,-c(1:3)]
totalarea <- sum(suit.area[1,]) #historical distribution 
small <- which(apply(suit.area, 2, sum, na.rm=T)/totalarea < 0.25) # establish insignificant species for removal

for(spp in spps[-small]){
  
  RefSuit <- get(paste("Suit.ref", spp, edatope, sep="."))
  outRange.base <- RefSuit==5
  RefSuit[RefSuit==4] <- 3 #added this based on email from Will May 18, 2021
  RefSuit[RefSuit==5] <- 4
  RefSuit[is.na(RefSuit)] <- 4
  table(RefSuit)
  
  for(proj.year in proj.years){
  # compile the GCM projections into a data frame
  ProjSuit <- data.frame(temp=rep(NA, length(RefSuit))) #initiate the data frame with a dummy column
  ChangeSuit <- data.frame(temp=rep(NA, length(RefSuit))) #initiate the data frame with a dummy column
  for(GCM in GCMs){
    temp <- get(paste("Suit", GCM, rcp, proj.year, spp, edatope, sep="."))
    # temp <- read.csv(paste(outdir, "/delete/Suit", studyarea, GCM, rcp, proj.year, spp, edatope, "csv", sep="."))[,1]
    temp[temp==4] <- 3 #added this based on email from Will May 18, 2021
    temp[temp==5] <- 4
    temp[is.na(temp)] <- 4
    ProjSuit <- cbind(ProjSuit,temp)
    ChangeSuit <- cbind(ChangeSuit,RefSuit-temp)
  }
  ProjSuit <- ProjSuit[,-1] #remove the dummy column
  ChangeSuit <- ChangeSuit[,-1] #remove the dummy column
  names(ProjSuit) <- GCMs
  names(ChangeSuit) <- GCMs
  
  # calculate ensemble mean suitability change. this isn't biased by missing suitabilties for exotic BGCs
  ChangeSuit.mean <- apply(ChangeSuit, 1, mean, na.rm=T)
  
  outRange <- outRange.base
  outRange[which(ChangeSuit.mean!=0)] <- FALSE
  ChangeSuit.mean[outRange==T] <- NA
  
  values(X) <- NA
  values(X)[land] <- ChangeSuit.mean
  # plot(X)
  writeRaster(X, paste(outdir, "/data/Spp.ChangeSuit", studyarea, spp, edatope, rcp, proj.year,"tif", sep="."),overwrite=TRUE, format="GTiff")
  
  ##=================================
  # binary appearance/disappearance
  Suit.ensemble <- as.matrix(ProjSuit)
  Suit.ensemble[Suit.ensemble==5] <- 4
  binary <- rep(0, length(RefSuit))
  binary[outRange.base==T] <- NA
  binary[outRange.base] <- apply(Suit.ensemble[outRange.base,], 1, function(x){return(if((sum(x<4, na.rm=T)/sum(!is.na(x)))>0) sum(x<4, na.rm=T)/sum(!is.na(x)) else NA)})
  binary[outRange.base==F] <- apply(Suit.ensemble[outRange.base==F,], 1, function(x){return(0-sum(x==4, na.rm=T)/sum(!is.na(x)))})
  values(X) <- NA
  values(X)[land] <- binary
  writeRaster(X, paste(outdir, "/data/Spp.binary", studyarea, spp, edatope, rcp, proj.year,"tif", sep="."),overwrite=TRUE, format="GTiff")

  # print(proj.year)
}
  
    print(spp)
}

print(edatope)
}

