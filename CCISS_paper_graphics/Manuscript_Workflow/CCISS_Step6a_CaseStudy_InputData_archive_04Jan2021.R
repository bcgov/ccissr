
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

P4S.latlon <- CRS("+proj=longlat +datum=WGS84")
P4S.albers <- CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs")

studyarea <- "Tahltan"

#==================================================
# Select the study area based on a bounding box around the focal BGC unit, and create the dem and ClimateWNA file from scratch
# Only need to do this for a new focal unit
#==================================================

#===============================================================================
# DEM and boundary
#===============================================================================

## boundary for the study area
bdy.aea <- readOGR("C:/GitHub/2019_CCISS/inputs/shapes/Tahltan/TahltanLUP_BoundaryV7_200702.shp")

bdy <- spTransform(bdy.aea, P4S.latlon) #reproject the gridpoints
plot(bdy)

# ####### create a polygon mask for the focal unit 
border=0.5
my_box = as(extent(buffer(bdy, border)), "SpatialPolygons")      		# convert extent box to shapefile (rectangle)
proj4string(my_box) = projection(bdy)				# assign spatial projection to extent object
focal.mask <- gDifference(my_box, bdy)
projection(focal.mask)  #verify latlong projection of the study area boundary
plot(focal.mask)
plot(my_box)
plot(bdy, add=T)

## create the study area
border <- 0.25 # buff
r <- raster("C:\\Users/CMAHONY/OneDrive - Government of BC\\SpatialData\\DEM\\namer_dem1.bil")
r <- crop(r, buffer(bdy, border))
plot(r)
plot(bdy, add=T)

writeRaster(r, filename=paste("inputs\\dem1_", studyarea, ".tif", sep=""), format="GTiff", overwrite=TRUE)
## create a data frame of the projected coordinates of the dem
r.pts <- rasterToPoints(r, spatial=T) #create a spatial points data frame of the non-NA values of the DEM
projection(r.pts) <- P4S.latlon 
str(r.pts)

### Attribute the points with BGC v10 linework
BGCv11 <- readOGR("C:\\Users\\CMAHONY/OneDrive - Government of BC\\SpatialData\\BGCv11\\BGCv11_WithLandcover.gdb")
projection(BGCv11) <- P4S.albers # tell R the projection of the shapefile (this is known via the shapefile metadata)
BGCv11 <- crop(BGCv11, buffer(bdy.aea, 40000))

pts.aea <- spTransform(r.pts, P4S.albers) # reproject to lat-long 
# plot(pts.aea, add=T, pch=16, cex=0.5)
BGC.pts <- over(pts.aea,BGCv11)  ##extract the BGCv11 attributes to the points (this takes a while, like 30ish minutes)
str(BGC.pts)


r.pts <- as.data.frame(r.pts) #projected coordinates of the dem

## create the climateNA input file
CNAinput <- data.frame(id1=1:dim(r.pts)[1], id2=as.character(BGC.pts$MAP_LABEL), lat=r.pts$y, lon=r.pts$x, el=r.pts$namer_dem1)
str(CNAinput)
write.csv(CNAinput,paste("inputs\\", studyarea, "Grid.csv", sep=""), row.names=FALSE)
#NOTE! FUTURE NORMALS SHOULD BE QUERIED FROM CLIMATEBC AS BATCH FILES FOR ALL 15 GCMS AND FOR ONLY ONE RCP AND NORMAL PERIOD. 


#===============================================================================
# Set analysis Parameters
#===============================================================================

source("./_CCISS_Packages.R") ## packages required
source("./_CCISS_Functions.R") ## common functions
source("./_CCISS_Parameters.R") ## settings used through all scripts

grid <- paste(studyarea, "Grid", sep="")

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
#rownames(importance(BGCmodel)) ### shows the variable used in the RFmodel
vars <- as.data.frame(BGCmodel$variable.importance)
vars <- row.names(vars)
# setwd('C:/GitHub/2019_CCISS')

Columns <- unique(c("PPT05", "PPT06", "PPT07", "PPT08", "PPT09", "PPT_at", 
                    "PPT_wt", "CMD07", "CMD", "MAT", "PPT_sm", "Tmin_wt", "Tmax_sm",
                    vars[!vars %in% c("PPT_MJ", "PPT_JAS", "PPT.dormant", "CMD.def", "CMDMax", "CMD.total")]))


#===============================================================================
# create a hillshade backdrop
#===============================================================================


dem <- raster(paste("inputs\\dem1_",studyarea,".tif", sep=""))

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
# alt <- crop(alt,X)
slope = terrain(alt, opt='slope')
aspect = terrain(alt, opt='aspect')
hill = hillShade(slope, aspect, 40, 270)
plot(hill, col=grey(0:100/100), legend=FALSE)

## subset the raster for focal cells only
focal.buffer <- gBuffer(focal.mask, width=-res(dem)[1])
focal.cells.buffer <- which(is.na(values(mask(dem,focal.buffer))))
focal.cells <- which(is.na(values(mask(dem,focal.mask))))

X <- dem
values(X) <- NA
values(X)[focal.cells] <- 1
plot(X)

plot(dem)
plot(bdy, add=T, border="blue")
plot(hill, col=grey(0:100/100), legend=FALSE)
plot(focal.mask, add=T, col=alpha("white", 0.5))

#===============================================================================
#===============================================================================
# BGC Projections
#===============================================================================
#===============================================================================

#===============================================================================
# BGC Projections for reference period
#===============================================================================

# setwd("C:\\Users/CMAHONY/OneDrive - Government of BC\\Projects\\2019_CCISS")

fplot=paste("inputs\\", grid, "_Normal_1961_1990MSY.csv", sep="")

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
BGC.pred.ref <- predict(BGCmodel, Y0)
write.csv(list(BGC.pred.ref$predictions), paste("outputs\\BGC.pred", grid, "ref.csv", sep="."), row.names = F)

# ## Write Climate file ######
# write.csv(Y0[,which(names(Y0)%in%VarList)], paste("inputs\\", grid, "_1961_1990_BioVars.csv", sep=""), row.names = F)

## calculate mean climate of study area for use in calculating change
Y.refmean <- apply(Y0, 2, FUN=mean, na.rm=T)[-c(1:5)]
write.csv(t(as.data.frame(Y.refmean)), paste("outputs\\clim.refMean",grid,"csv", sep="."), row.names = F)

change <- data.frame("GCM"="obs", "rcp"="obs", "proj.year"=1975, as.data.frame(t(rep(0, length(Y.refmean)))))
names(change)[-c(1:3)] <- names(Y.refmean)

dem <- raster(paste("inputs\\dem1_",studyarea,".tif", sep=""))
instudyarea <- unlist(cellFromPolygon(dem, bdy))
land <- which(!is.na(values(dem)[instudyarea]))
X <- dem
values(X) <- NA
values(X)[land] <- BGC.pred.ref$predictions
writeRaster(X, paste("Deliverables\\Novelty", analogPool, RCP, proj.year,"tif", sep="."),overwrite=TRUE, format="GTiff")

plot(X)

#===============================================================================
# BGC Projections for recent period
#===============================================================================

# setwd("C:\\Users/CMAHONY/OneDrive - Government of BC\\Projects\\2019_CCISS")
hist.years <- c(1995, 2005, 2015)
hist.periods <- c("1991_2000", "2001_2010", "2011_2019")

# read in the data for the decades
for(hist.year in hist.years){
  hist.period <- hist.periods[which(hist.years==hist.year)]
  fplot=paste("inputs\\", grid, "_Decade_", hist.period, "MSY.csv", sep="")

  Y0 <- fread(fplot, stringsAsFactors = FALSE, data.table = FALSE) #fread is faster than read.csv

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
str(Y.2005)

## take the weighted mean of the decades
Y.hist <- (Y.1995*10+Y.2005*10+Y.2015*9)/29

##Predict future subzones######
BGC.pred.hist <- predict(BGCmodel, Y.hist)
write.csv(list(BGC.pred.hist$predictions), paste("outputs\\BGC.pred", grid, "2005.csv", sep="."), row.names = F)

## calculate climate change
Y.mean <- apply(Y.hist, 2, FUN=mean, na.rm=T)[-c(1:5)]
change.temp <- Y.mean - Y.refmean
change <- rbind(change, data.frame("GCM"="obs", "rcp"="obs", "proj.year"=2005, as.data.frame(t(change.temp))))

#===============================================================================
# BGC Projections for future periods
#===============================================================================
rcp="rcp45"

# for(rcp in rcps){
  for(proj.year in proj.years){
    fplot=paste("inputs\\", grid, "_", rcp, "_", proj.year, "_MSY.csv", sep="")
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
      BGC.pred.proj <- predict(BGCmodel, temp)
      BGC.pred.proj <- list(BGC.pred.proj$predictions)
      write.csv(BGC.pred.proj, paste("outputs\\BGC.pred",grid, GCM, rcp, proj.year,"csv", sep="."), row.names = F)

      ## calculate climate change
      Y.mean <- apply(temp[,which(names(temp)%in%names(Y.refmean))], 2, FUN=mean, na.rm=T)
      
      change.temp <- Y.mean - Y.refmean
      change <- rbind(change, data.frame("GCM"=GCM, "rcp"=rcp, "proj.year"=proj.year, as.data.frame(t(change.temp))))
      
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

write.csv(change, paste("outputs\\clim.meanChange",grid,"csv", sep="."), row.names = F)

#===============================================================================
# Ensemble Mean BGC Projections for future periods
#===============================================================================

# determine vote winner BGC and ensemble agreement (WARNING: takes about 1 minute per rcp/proj.year)
# for(rcp in rcps){
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
    write.csv(get(paste("BGC.pred.ensemble", rcp, proj.year, sep=".")), paste("outputs\\BGC.pred",grid, "ensemble", rcp, proj.year,"csv", sep="."), row.names = F)
    # write.csv(agreement, paste("outputs\\BGC.pred",grid, "agreement", rcp, proj.year,"csv", sep="."), row.names = F)
    print(proj.year)
  }
  # print(rcp)
# }


#===============================================================================
# Make and export summary tables of BGC units for each future
#===============================================================================


BGC.pred.ref <- read.csv(paste("outputs\\BGC.pred", grid, "ref.csv", sep="."))[,1]
PredSum <- data.frame("GCM"="obs", "rcp"="obs", "proj.year"=1975, as.data.frame(table(BGC.pred.ref)))
names(PredSum)[4] <- "BGC.pred"

# Historical BGC
hist.years=2005
hist.year=2005
for(hist.year in hist.years){
  BGC.pred <- read.csv(paste("outputs\\BGC.pred", grid, "2005.csv", sep="."))[,1]
  assign(paste("BGC.pred", hist.year, sep="."), BGC.pred) #bgc projection
  PredSum <- rbind(PredSum, data.frame("GCM"="obs", "rcp"="obs", "proj.year"=2005, as.data.frame(table(BGC.pred))))
  print(hist.year)
}

# Future BGC
# for(rcp in rcps){
for(proj.year in proj.years){
  for(GCM in c("ensemble", GCMs)){
    BGC.pred <- as.character(read.csv(paste("outputs\\BGC.pred",grid, GCM, rcp, proj.year,"csv", sep="."), header = T)[,1])
    assign(paste("BGC.pred", GCM, rcp, proj.year, sep="."), BGC.pred) #bgc projection
    PredSum <- rbind(PredSum, data.frame("GCM"=GCM, "rcp"=rcp, "proj.year"=proj.year, as.data.frame(table(BGC.pred))))
    # print(GCM)
  }
  print(proj.year)
}
# print(rcp)
# }

# write out summary of BGC units for each future. 
PredSum.wide <- reshape(PredSum, idvar = c(names(PredSum)[1:3]), timevar = "BGC.pred", direction = "wide")
write.csv(PredSum.wide, paste("outputs\\PredSum.BGC",grid,"csv", sep="."), row.names = F)

#===============================================================================
#===============================================================================
# Species Feasibility Projections
#===============================================================================
#===============================================================================

# #===============================================================================
# # generate the vector of mapped BGCs
# #===============================================================================
# 
# points <- fread(paste("./inputs/", grid,".csv", sep=""))
# BGC <- points$ID2
# BGC <- gsub(" ","",BGC)  
# sort(table(BGC))

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


#===============================================================================
# find the species suitability each projection/edatope/species combination
#===============================================================================

# Import suitability tables
S1 <- treesuit
S1 <- unique(S1)[,1:4]
dim(S1)
S1 <- unique(S1)
dim(S1)

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
    BGC.pred <- as.character(BGC.pred.ref) # get the BGC prediction
    
    # get the suitability for the selected species associated with each site series
    suit <- S1$ESuit[which(S1$Spp==spp)][match(as.vector(unlist(SiteLookup[which(names(SiteLookup)==edatope)])), as.character(S1$SS_NoSpace[which(S1$Spp==spp)]))]
    Suit.ref <- suit[match(BGC.pred, SiteLookup$BGC)]
    Suit.ref[is.na(Suit.ref)] <- 5 #set the NA values to suitability 5
    Suit.ref[Suit.ref==4] <- 5 #set 4 to suitability 5
    write.csv(Suit.ref, paste("outputs\\Suit.ref", grid, spp, edatope, "csv", sep="."), row.names = F)
    
    for(hist.year in hist.years){
      BGC.pred <- as.character(get(paste("BGC.pred", hist.year, sep=".")))
      ## identify cells with no suitability interpretations
      bgc.exotic <- (1:length(BGC.pred))[-which(BGC.pred%in%unique(BGC))]
      bgc.exotic.noSuit <- bgc.exotic[-which(BGC.pred[bgc.exotic]%in%unique(S1$BGC))]
      
      suit <- S1$ESuit[which(S1$Spp==spp)][match(as.vector(unlist(SiteLookup[which(names(SiteLookup)==edatope)])), S1$SS_NoSpace[which(S1$Spp==spp)])]
      temp <- suit[match(BGC.pred, SiteLookup$BGC)]
      temp[is.na(temp)] <- 5 #set the NA values to suitability 5
      temp[temp==4] <- 5 #set 4 to suitability 5
      temp[bgc.exotic.noSuit] <- NA # set cells with no suitabilty interpretatoin to NA
      assign(paste("Suit", hist.year, sep="."), temp)
      write.csv(temp, paste("outputs\\Suit", grid, hist.year, spp, edatope, "csv", sep="."), row.names = F)
      # print(hist.year)
    }
    
    # get the suitability for future periods, for each projection/edatope/species combination
    for(GCM in c(GCMs, "ensemble")){
      # for(rcp in rcps){
        for(proj.year in proj.years){
          # get the BGC projection and sub in the crosswalk between the modeled units and the table units
          BGC.pred <- as.character(read.csv(paste("outputs\\BGC.pred",grid, GCM, rcp, proj.year,"csv", sep="."), header = T)[,1])
          # BGC.pred[which(BGC.pred%in%Crosswalk$Modeled)] <- as.character(Crosswalk$Tables[match(BGC.pred[which(BGC.pred%in%Crosswalk$Modeled)], Crosswalk$Modeled)]) # XXX THIS IS NOT CORRECT. NEED TO FIGURE OUT HOW TO INCORPORATE THE CROSSWALK TABLE PROPERLY. sub in the crosswalk between the modeled units and the table units

          ## identify cells with no suitability interpretations
          bgc.exotic <- (1:length(BGC.pred))[-which(BGC.pred%in%unique(BGC))]
          bgc.exotic.noSuit <- bgc.exotic[-which(BGC.pred[bgc.exotic]%in%unique(S1$BGC))]

          # get the suitability for the selected species associated with each site series
          suit <- S1$ESuit[which(S1$Spp==spp)][match(as.vector(unlist(SiteLookup[which(names(SiteLookup)==edatope)])), S1$SS_NoSpace[which(S1$Spp==spp)])]
          temp <- suit[match(BGC.pred, SiteLookup$BGC)]
          temp[is.na(temp)] <- 5 #set the NA values to suitability 5 (weights unsuitable a bit more heavily than suitable classes during averaging)
          temp[temp==4] <- 5 #set 4 to suitability 5
          temp[bgc.exotic.noSuit] <- NA # set cells with no suitabilty interpretation to NA
          assign(paste("Suit", GCM, rcp, proj.year, sep="."), temp)
          write.csv(temp, paste("outputs\\Suit", grid, GCM, rcp, proj.year, spp, edatope, "csv", sep="."), row.names = F)
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
points <- read.csv(paste("inputs\\",grid,".csv", sep=""))
BGC <- points$id2
BGC <- gsub(" ","",BGC)
exclude <- which(BGC%in%BGCs_notin_THLB$BGC[which(BGCs_notin_THLB$Exlude=="x")])


for(edatope in edatopes){
  
  #initiate tables to store summary values
  PredSum.suit <- data.frame(PredSum.wide[1:3], as.data.frame(matrix(rep(NA, length(spps)*dim(PredSum.wide)[1]), dim(PredSum.wide)[1])))
  names(PredSum.suit)[-c(1:3)] <- spps
  PredSum.spp <- PredSum.suit
  PredSum.suit.home <- PredSum.suit #home is for counting cells within historical range. 
  PredSum.spp.home <- PredSum.suit
  
  for(spp in spps){
    
    #reference period suitabilities
    Suit.ref <- read.csv(paste("outputs\\Suit.ref", grid, spp, edatope, "csv", sep="."))[-exclude,1]
    # Suit.ref <- read.csv(paste("outputs\\Suit.ref", grid, spp, edatope, "csv", sep="."))[,1]
    Suit.ref[Suit.ref==5] <- NA
    outRange.ref <- is.na(Suit.ref)
    Suit.ref[is.na(Suit.ref)] <- 5
    Suit.ref <- 1-(Suit.ref-1)/4

    row <- 1
    col <- which(names(PredSum.suit)==spp)
    PredSum.suit[row,col] <- round(sum(Suit.ref))
    PredSum.spp[row,col] <- round(sum(Suit.ref>0))
    PredSum.suit.home[row,col] <- round(sum(Suit.ref[outRange.ref==F]))
    PredSum.spp.home[row,col] <- round(sum((Suit.ref>0)[outRange.ref==F]))
    
    for(hist.year in hist.years){
      Suit.proj <- read.csv(paste("outputs\\Suit", grid, hist.year, spp, edatope, "csv", sep="."))[-exclude,1]
      Suit.proj[is.na(Suit.proj)] <- 5
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
          Suit.proj <- read.csv(paste("outputs\\Suit", grid, GCM, rcp, proj.year, spp, edatope, "csv", sep="."))[-exclude,1]
          Suit.proj[is.na(Suit.proj)] <- 5
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

  write.csv(PredSum.suit, paste("outputs\\PredSum.suit",grid, edatope,"csv", sep="."), row.names = F)
  write.csv(PredSum.spp, paste("outputs\\PredSum.spp",grid, edatope,"csv", sep="."), row.names = F)
  write.csv(PredSum.suit.home, paste("outputs\\PredSum.suit.home",grid, edatope,"csv", sep="."), row.names = F)
  write.csv(PredSum.spp.home, paste("outputs\\PredSum.spp.home",grid, edatope,"csv", sep="."), row.names = F)
  
  
    print(edatope)
}

