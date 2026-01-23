
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


focal.unit <- "IDFdk3"


#==================================================
# Select the study area based on a bounding box around the focal BGC unit, and create the dem and ClimateWNA file from scratch
# Only need to do this for a new focal unit
#==================================================

setwd("C:\\Colin\\Projects\\2019_CCISS")

grid <- "BC2kmGrid"

## get the vector of BGC units
points <- read.csv(paste("InputData\\",grid,".csv", sep=""))
BGC <- as.character(points$ID2)
BGC <- gsub(" ","",BGC)

## find the bounding box for the focal unit
focal.lat <- range(points$lat[which(BGC==focal.unit)])
focal.lon <- range(points$lon[which(BGC==focal.unit)])

## create the study area
border <- 0.25 # buff
r <- raster("C:\\Colin\\SpatialData\\DEM\\namer_dem1.bil")
plot(r, xlim=focal.lon+c(-border, border), ylim=focal.lat+c(-border/2, border/2))
rect(focal.lon[1],focal.lat[1],focal.lon[2],focal.lat[2])

r.crop <- crop(r, c(focal.lon+c(-border, border), focal.lat+c(-border/2, border/2)))
# r.crop <- aggregate(r.crop, fact=2)
writeRaster(r.crop, filename=paste("InputData\\dem1_", focal.unit, ".tif", sep=""), format="GTiff", overwrite=TRUE)
## create a data frame of the projected coordinates of the dem
r.pts <- rasterToPoints(r.crop, spatial=T) #create a spatial points data frame of the non-NA values of the DEM
r.pts <- as.data.frame(r.pts) #projected coordinates of the dem

## create the climateNA input file
CNAinput <- data.frame(id1=1:dim(r.pts)[1], id2=rep(NA, dim(r.pts)[1]), lat=r.pts$y, lon=r.pts$x, el=r.pts$namer_dem1)
str(CNAinput)
write.csv(CNAinput,paste("InputData\\", focal.unit, "Grid.csv", sep=""), row.names=FALSE)
#NOTE! FUTURE NORMALS SHOULD BE QUERIED FROM CLIMATEBC AS BATCH FILES FOR ALL 15 GCMS AND FOR ONLY ONE RCP AND NORMAL PERIOD. 

#===============================================================================
# Set analysis Parameters
#===============================================================================

# setwd("C:\\Colin\\Projects\\2019_CCISS")

grid <- paste(focal.unit, "Grid", sep="")

GCMs <-  c("ACCESS1-0","CanESM2","CCSM4","CESM1-CAM5","CNRM-CM5","CSIRO-Mk3-6-0", "GFDL-CM3","GISS-E2R", "HadGEM2-ES", "INM-CM4", "IPSL-CM5A-MR", "MIROC-ESM", "MIROC5", "MPI-ESM-LR","MRI-CGCM3")    
# GCMs <-  GCMs[1:2] #XXX for testing purposes

##Specify variables in model
model = "5.1"
VarList = c("AHM", "bFFP","CMD.total","DD5_sp","EMT","Eref_sm","EXT","FFP","MCMT","MSP",
            "PPT_JAS","PPT_MJ","PPT06","SHM","TD","Tmax_sp","Tmin_at","Tmin_sm","Tmin_wt",
            "PAS","CMD.def","CMDMax","eFFP","Eref09","MAT","PPT07","Tmin_sp")

###Load random forest model
fname="Rcode\\BGCv11_AB_USA_16VAR_SubZone_RFmodel.Rdata"
#fname = (file.choose())
load(fname)
VarList = rownames(importance(BGCmodel))

GCMs <-  c("ACCESS1-0","CanESM2","CCSM4","CESM1-CAM5","CNRM-CM5","CSIRO-Mk3-6-0", "GFDL-CM3","GISS-E2R", "HadGEM2-ES", "INM-CM4", "IPSL-CM5A-MR", "MIROC-ESM", "MIROC5", "MPI-ESM-LR","MRI-CGCM3")
rcps <- c("rcp45", "rcp85")
proj.years <- c(2025, 2055, 2085)
hist.years <- c(1995, 2004, 2005, 2009, 2014, 2017)
edatopes<- c("B2", "C4", "D6")
spps.lookup <- read.csv("InputData\\Tree speciesand codes_2.0.csv")
edatope.name <- c("Subxeric-poor", "Mesic-medium", "Hygric-rich")
BGCcolors <- read.csv("C:\\Colin\\Projects\\2019_CCISS\\InputData\\BGCzone_Colorscheme.csv")

#===============================================================================
# DEM and boundary
#===============================================================================
dem <- raster(paste("InputData\\dem1_",focal.unit,".tif", sep=""))

## read in boundary shapefile (created separately)
P4S.AEA <- CRS ("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs") # standard albers projection for BC gov't data
focal.bdy <- readShapePoly(paste("InputData\\", focal.unit, "_BEC11.shp", sep=""))
projection(focal.bdy) <- P4S.AEA
#reproject shapefiles to dem projection
P4S.latlon <- CRS("+proj=longlat +datum=WGS84")
focal.bdy <- spTransform(focal.bdy, P4S.latlon) # reproject to lat-long 

# ####### create a polygon mask for the focal unit 
border=0.5
my_box = as(extent(c(focal.lon+c(-border, border), focal.lat+c(-border/2, border/2))), "SpatialPolygons")      		# convert extent box to shapefile (rectangle)
proj4string(my_box) = projection(focal.bdy)				# assign spatial projection to extent object
focal.mask <- gDifference(my_box, focal.bdy)
projection(focal.mask)  #verify latlong projection of the study area boundary

X <- dem

## create a hillshade backdrop
mapsheets <- c("92I|92J|92O|92P|93A|93B")
setwd("C:\\Colin\\Projects\\SpatialData\\CDEM")
files <- list.files(pattern="*.tif")
files <- files[grep(mapsheets, files)]
for(i in 1:length(files)){
  temp <- aggregate(raster(files[i]), fact=8) #reduce resolution
  if(i==1) alt <- temp else alt <- mosaic(alt, temp, fun=mean)
  print(i)
}
setwd("C:\\Colin\\Projects\\2019_CCISS")
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
plot(focal.bdy, add=T, border="blue")
plot(hill, col=grey(0:100/100), legend=FALSE)
plot(focal.mask, add=T, col=alpha("white", 0.5))


#===============================================================================
# BGC Projections for reference period
#===============================================================================

# setwd("C:\\Colin\\Projects\\2019_CCISS")

fplot=paste("InputData\\", grid, "_Normal_1961_1990MSY.csv", sep="")

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
write.csv(BGC.pred.ref, paste("OutputData\\BGC.pred.ref", grid, "csv", sep="."), row.names = F)

## Write Climate file ######
write.csv(Y0[,which(names(Y0)%in%VarList)], paste("InputData\\", grid, "_1961_1990_BioVars.csv", sep=""), row.names = F)

#===============================================================================
# BGC Projections for historical decades
#===============================================================================

# setwd("C:\\Colin\\Projects\\2019_CCISS")
hist.years <- c(1995, 2005)
hist.periods <- c("1991_2000", "2001_2010")

for(hist.year in hist.years){
  hist.period <- hist.periods[which(hist.years==hist.year)]
  fplot=paste("InputData\\", grid, "_Decade_", hist.period, "MSY.csv", sep="")

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
  assign(paste("BGC.pred",hist.year,sep="."), predict(BGCmodel, Y0))
  write.csv(get(paste("BGC.pred",hist.year,sep=".")), paste("OutputData\\BGC.pred", grid,hist.year,"csv", sep="."), row.names = F)

  ## Write Climate file ######
  write.csv(Y0[,which(names(Y0)%in%VarList)], paste("InputData\\", grid, "_",hist.year,"_BioVars.csv", sep=""), row.names = F)

  print(hist.year)
}

#===============================================================================
# BGC Projections for last decade
#===============================================================================

# setwd("C:\\Colin\\Projects\\2019_CCISS")
  fplot=paste("InputData\\", grid, "_2011-2017MSYT.csv", sep="")

  Y0 <- fread(fplot, stringsAsFactors = FALSE, data.table = FALSE) #fread is faster than read.csv
  # Y0 <- Y0[!is.na(Y0[,2]),]
str(Y0)

  #####generate some additional variables
  Y0$PPT_MJ <- Y0$PPT05 + Y0$PPT06 # MaY/June precip
  Y0$PPT_JAS <- Y0$PPT07 + Y0$PPT08 + Y0$PPT09 # July/Aug/Sept precip
  Y0$PPT.dormant <- Y0$PPT_at + Y0$PPT_wt # for calculating spring deficit
  Y0$CMD.def <- 500 - (Y0$PPT.dormant)# start of growing season deficit original value was 400 but 500 seems better
  Y0$CMD.def [Y0$CMD.def < 0] <- 0 #negative values set to zero = no deficit
  Y0$CMDMax <- Y0$CMD07
  Y0$CMD.total <- Y0$CMD.def + Y0$CMD

  # Extract the year 2017
  Y2017 <- Y0[which(Y0$Year==2017),]
str(Y2017)

  # Calculate mean of 2011-2017 period
  Y1 <- aggregate(Y0, by=list(Y0$id1), FUN=mean, na.rm=T)[,-1]
  Y1 <- Y1[match(Y2017$id1, Y1$id1),]

  ##Predict BGC units######
  BGC.pred.2017 <- predict(BGCmodel, Y2017)
  BGC.pred.2014 <- predict(BGCmodel, Y1)
  write.csv(BGC.pred.2017, paste("OutputData\\BGC.pred", grid,"2017.csv", sep="."), row.names = F)
  write.csv(BGC.pred.2014, paste("OutputData\\BGC.pred", grid,"2014.csv", sep="."), row.names = F)

  ## Write Climate file ######
  write.csv(Y2017[,which(names(Y2017)%in%VarList)], paste("InputData\\", grid, "_2017_BioVars.csv", sep=""), row.names = F)
  write.csv(Y1[,which(names(Y1)%in%VarList)], paste("InputData\\", grid, "_2014_BioVars.csv", sep=""), row.names = F)

# ===============================================================================
# BGC Projections for other historical normals
# ===============================================================================

# setwd("C:\\Colin\\Projects\\2019_CCISS")
hist.years <- c(1995, 2005)
hist.periods <- c("1991_2000", "2001_2010")

mean(c(1991,2017))
mean(c(2001,2017))

# read in the data for the 1990s and 2000s
for(hist.year in hist.years){
  hist.period <- hist.periods[which(hist.years==hist.year)]
  fplot=paste("InputData\\", grid, "_Decade_", hist.period, "MSY.csv", sep="")

  Y0 <- fread(fplot, stringsAsFactors = FALSE, data.table = FALSE) #fread is faster than read.csv

  #####generate some additional variables
  Y0$PPT_MJ <- Y0$PPT05 + Y0$PPT06 # MaY/June precip
  Y0$PPT_JAS <- Y0$PPT07 + Y0$PPT08 + Y0$PPT09 # July/Aug/Sept precip
  Y0$PPT.dormant <- Y0$PPT_at + Y0$PPT_wt # for calculating spring deficit
  Y0$CMD.def <- 500 - (Y0$PPT.dormant)# start of growing season deficit original value was 400 but 500 seems better
  Y0$CMD.def [Y0$CMD.def < 0] <- 0 #negative values set to zero = no deficit
  Y0$CMDMax <- Y0$CMD07
  Y0$CMD.total <- Y0$CMD.def + Y0$CMD
  assign(paste("Y",hist.year,sep="."), Y0[,which(names(Y0)%in%VarList)])
  print(hist.year)
}
str(Y.2005)

# 2011-2017 period already exported in previous phase, so read that in
Y.2014 <- read.csv(paste("InputData\\", grid, "_2014_BioVars.csv", sep=""))
str(Y.2014)

Y.2009 <- (Y.2005*10+Y.2014*7)/17
Y.2004 <- (Y.1995*10+Y.2005*10+Y.2014*7)/27
str(Y.2009)
str(Y.2004)

hist.years <- c(2004, 2009)
for(hist.year in hist.years){

  ##Predict future subzones######
  assign(paste("BGC.pred",hist.year,sep="."), predict(BGCmodel, get(paste("Y",hist.year,sep="."))))
  write.csv(get(paste("BGC.pred",hist.year,sep=".")), paste("OutputData\\BGC.pred", grid,hist.year,"csv", sep="."), row.names = F)

  ## Write Climate file ######
  write.csv(get(paste("Y",hist.year,sep=".")), paste("InputData\\", grid, "_",hist.year,"_BioVars.csv", sep=""), row.names = F)

  print(hist.year)
}

#===============================================================================
# BGC Projections for future periods
#===============================================================================
rcp="rcp45"

# for(rcp in rcps){
  for(proj.year in proj.years){
    fplot=paste("InputData\\", grid, "_", rcp, "_", proj.year, "_MSYT.csv", sep="")
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
    Ystr <- strsplit(Y0[,1], "_")
    Y4 <- matrix(unlist(Ystr), ncol=3, byrow=TRUE)
    
    GCMs <- unique(Y4[,1])
    for(GCM in GCMs){
      temp <- Y0[which(Y4[,1]==GCM & Y4[,2]==rcp & Y4[,3]==proj.year),]
      assign(paste("BGC.pred", GCM, rcp, proj.year, sep="."), predict(BGCmodel, temp))
      write.csv(get(paste("BGC.pred", GCM, rcp, proj.year, sep=".")), paste("OutputData\\BGC.pred",grid, GCM, rcp, proj.year,".csv", sep=""), row.names = F)
      print(GCM)
    }
    print(proj.year)
  }
  # print(rcp)
# }

#===============================================================================
# Import BGC projections for each period
#===============================================================================
hist.years <- c(1995, 2004, 2005, 2009, 2014, 2017)
hist.year.name <- c("1991-2000", "1991-2017", "2001-2010", "2001-2017","2011-2017", "2017")

## reference period BGC
BGC.pred.ref <- as.character(read.csv(paste("OutputData\\BGC.pred.ref", grid, "csv", sep="."))[,1])
zone.pred.ref <- rep(NA, length(BGC.pred.ref))
for(i in BGCcolors$zone){ zone.pred.ref[grep(i,BGC.pred.ref)] <- i }

# Historical BGC
for(hist.year in hist.years){
  BGC.pred <- as.character(read.csv(paste("OutputData\\BGC.pred", grid,hist.year,"csv", sep="."))[,1])
  assign(paste("BGC.pred",hist.year,sep="."), BGC.pred)
  print(hist.year)
}

# Future BGC
PredSum <- data.frame()
for(rcp in rcps){
  for(proj.year in proj.years){
    for(GCM in GCMs){
      BGC.pred <- as.character(read.csv(paste("OutputData\\BGC.pred",grid, GCM, rcp, proj.year,".csv", sep=""))[,1])
      assign(paste("BGC.pred", GCM, rcp, proj.year, sep="."), BGC.pred)
      # print(GCM)
      PredSum <- rbind(PredSum, as.data.frame(table(BGC.pred)))
    }
    print(proj.year)
  }
  print(rcp)
}
PredSum <- aggregate(PredSum$Freq, by=list(PredSum$BGC.pred), FUN=sum)
names(PredSum) <- c("BGC", "count")


#===============================================================================
# import will's lookup table for the site series associated with the selected edatope in each BGC unit
#===============================================================================
SiteSeries_Use <-read.csv(paste("InputData/","SiteSeries_Use_5",".csv",sep=""),stringsAsFactors=FALSE,na.strings=".")
str(SiteSeries_Use)

SiteLookup <- data.frame(MergedBGC=unique(SiteSeries_Use$MergedBGC))
edatopes<- c("B2", "C4", "D6")
for(edatope in edatopes){
  # SiteLookup <- cbind(SiteLookup, SiteSeries_Use$SS_NoSpace[match(SiteLookup[,1], SiteSeries_Use$MergedBGC[which(SiteSeries_Use$Use==edatope)])])
  SiteLookup <- cbind(SiteLookup, SiteSeries_Use$SS_NoSpace[which(SiteSeries_Use$Edatopic==edatope)])
  names(SiteLookup)[which(edatopes==edatope)+1] <- edatope
}
str(SiteLookup)

#===============================================================================
# find the species suitability each projection/edatope/species combination
#===============================================================================

# Import suitability tables
treesuit="TreeSpp_ESuit_v11_17"
S1 <- read.csv(paste("InputData/",treesuit,".csv",sep=""),stringsAsFactors=F,na.strings=".")
S1 <- unique(S1)

## EDA: are there suitabilities for all projected units? 
NoSuit <- PredSum[-which(PredSum$BGC%in%S1$BGC),]
NoSuit[rev(order(NoSuit$count)),]

## EDA: Which site series are missing suitabilities? 
for(edatope in edatopes) assign(paste("NoSuit", edatope, sep="."), SiteLookup[-which(SiteLookup[,which(names(SiteLookup)==edatope)]%in%S1$Unit),which(names(SiteLookup)==edatope)])
for(edatope in edatopes) print(get(paste("NoSuit", edatope, sep=".")))

## EDA: are there any units missing from the SiteSeries_Use table? 
BGClist <- unique(S1$BGC)
BGClist[-which(BGClist%in%SiteLookup$MergedBGC)]

# species list
spps <- unique(S1$Spp)
spps <- spps[-which(spps=="X")]

for(spp in spps){
  for(edatope in edatopes){
    # get the suitability for the reference period predicted BGC.
    BGC.pred <- as.character(BGC.pred.ref) # get the BGC prediction
    # BGC.pred[which(BGC.pred%in%Crosswalk$Modeled)] <- as.character(Crosswalk$Tables[match(BGC.pred[which(BGC.pred%in%Crosswalk$Modeled)], Crosswalk$Modeled)]) # XXX THIS IS NOT CORRECT. NEED TO FIGURE OUT HOW TO INCORPORATE THE CROSSWALK TABLE PROPERLY. sub in the crosswalk between the modeled units and the table units
    # get the suitability for the selected species associated with each site series
    suit <- S1$ESuit[which(S1$Spp==spp)][match(as.vector(unlist(SiteLookup[which(names(SiteLookup)==edatope)])), S1$Unit[which(S1$Spp==spp)])]
    Suit.ref <- suit[match(BGC.pred, SiteLookup$MergedBGC)]
    Suit.ref[is.na(Suit.ref)] <- 5 #set the NA values to suitability 5
    Suit.ref[Suit.ref==4] <- 5 #set 4 to suitability 5
    write.csv(Suit.ref, paste("OutputData\\Suit.ref", grid, spp, edatope, "csv", sep="."), row.names = F)
    
    for(hist.year in hist.years){
      BGC.pred <- as.character(get(paste("BGC.pred", hist.year, sep=".")))
      ## identify cells with no suitability interpretations
      bgc.exotic <- (1:length(BGC.pred))[-which(BGC.pred%in%unique(BGC))]
      bgc.exotic.noSuit <- bgc.exotic[-which(BGC.pred[bgc.exotic]%in%unique(S1$BGC))]
      
      suit <- S1$ESuit[which(S1$Spp==spp)][match(as.vector(unlist(SiteLookup[which(names(SiteLookup)==edatope)])), S1$Unit[which(S1$Spp==spp)])]
      temp <- suit[match(BGC.pred, SiteLookup$MergedBGC)]
      temp[is.na(temp)] <- 5 #set the NA values to suitability 5
      temp[temp==4] <- 5 #set 4 to suitability 5
      temp[bgc.exotic.noSuit] <- NA # set cells with no suitabilty interpretatoin to NA
      assign(paste("Suit", hist.year, sep="."), temp)
      write.csv(temp, paste("OutputData\\Suit", grid, hist.year, spp, edatope, "csv", sep="."), row.names = F)
    }
    
    # get the suitability for future periods, for each projection/edatope/species combination
    rcp=rcps[1]  
    # for(rcp in rcps){
      for(proj.year in proj.years[-3]){
        for(GCM in GCMs){
          # get the BGC projection and sub in the crosswalk between the modeled units and the table units
          BGC.pred <- as.character(get(paste("BGC.pred", GCM, rcp, proj.year, sep=".")))
          # BGC.pred[which(BGC.pred%in%Crosswalk$Modeled)] <- as.character(Crosswalk$Tables[match(BGC.pred[which(BGC.pred%in%Crosswalk$Modeled)], Crosswalk$Modeled)]) # XXX THIS IS NOT CORRECT. NEED TO FIGURE OUT HOW TO INCORPORATE THE CROSSWALK TABLE PROPERLY. sub in the crosswalk between the modeled units and the table units
          
          ## identify cells with no suitability interpretations
          bgc.exotic <- (1:length(BGC.pred))[-which(BGC.pred%in%unique(BGC))]
          bgc.exotic.noSuit <- bgc.exotic[-which(BGC.pred[bgc.exotic]%in%unique(S1$BGC))]
          
          # get the suitability for the selected species associated with each site series
          suit <- S1$ESuit[which(S1$Spp==spp)][match(as.vector(unlist(SiteLookup[which(names(SiteLookup)==edatope)])), S1$Unit[which(S1$Spp==spp)])]
          temp <- suit[match(BGC.pred, SiteLookup$MergedBGC)]
          temp[is.na(temp)] <- 5 #set the NA values to suitability 5 (weights unsuitable a bit more heavily than suitable classes during averaging)
          temp[temp==4] <- 5 #set 4 to suitability 5
          temp[bgc.exotic.noSuit] <- NA # set cells with no suitabilty interpretatoin to NA
          assign(paste("Suit", GCM, rcp, proj.year, sep="."), temp)
          write.csv(temp, paste("OutputData\\Suit", grid, GCM, rcp, proj.year, spp, edatope, "csv", sep="."), row.names = F)
          # print(GCM)
        }
        # print(proj.year)
      }
      # print(rcp)
    # }
    
    # print(edatope)
  }
  print(paste(spp, " (", round(which(spps==spp)/length(spps)*100, 0), "%)", sep=""))
}

# ============================================================
# ============================================================
# RESULTS
# ============================================================
# ============================================================

# ============================================================
# BGC projections
# ============================================================

# determine vote winner BGC and ensemble agreement (WARNING: takes about 2 minutes per rcp/proj.year)
rcp=rcps[1]
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

############## compile the color scheme
BGC.pred.ref.WNA <- read.csv(paste("OutputData\\BGC.pred.ref.WNA2.csv", sep="."))[,1]
rcp="rcp45"
proj.year=2055
BGC.pred.WNA <- read.csv(paste("OutputData\\BGC.pred.WNA2", rcp, proj.year, "csv", sep="."))[,1]

zone.pred.ref.WNA <- gsub("[:a-z:]","",BGC.pred.ref.WNA) 
zone.pred.ref.WNA <- gsub("[:1-9:]","",zone.pred.ref.WNA) 

zone.pred.WNA <- gsub("[:a-z:]","",BGC.pred.WNA) 
zone.pred.WNA <- gsub("[:1-9:]","",zone.pred.WNA) 

zones <- sort(unique(c(zone.pred.ref.WNA, zone.pred.WNA)))
zones <- factor(zones, levels=zones)

#BGC zones
BGCcolors <- read.csv("C:\\Colin\\Projects\\2019_CCISS\\InputData\\BGCzone_Colorscheme.csv")
colors = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)][-1]

ColScheme <- sample(colors, dim(PredSum)[1])
units <- sort(as.character(PredSum[,1]))

ColScheme.zone <- rep(NA, length(zones))
ColScheme.zone[which(zones%in%BGCcolors$zone)] <- as.character(BGCcolors$HEX[match(zones[which(zones%in%BGCcolors$zone)], BGCcolors$zone)])
set.seed(8)
ColScheme.zone[-which(zones%in%BGCcolors$zone)] <- sample(colors, length(zones[-which(zones%in%BGCcolors$zone)]))
ColScheme.zone <- factor(ColScheme.zone, levels=ColScheme.zone)

################# 15-panel maps of all models for one time period
X <- dem
GCM=GCMs[1]
rcp=rcps[1]
proj.year=proj.years[2]

png(filename=paste("Results\\CCISS", focal.unit, "BGCprojections", rcp, proj.year,"png",sep="."), type="cairo", units="in", width=10.5, height=6.5, pointsize=10, res=600)
par(mar=c(0.1,0.1, 0.1,0.1), mgp=c(2,0.25,0), mfrow=c(3,5))

for(GCM in GCMs){
  
  BGC.pred <- get(paste("BGC.pred", GCM, rcp, proj.year, sep="."))
  zone.pred <- rep(NA, length(BGC.pred))
  for(i in zones){ zone.pred[grep(i,BGC.pred)] <- i }
  bgcs <- table(BGC.pred[focal.cells])
  bgcs <- bgcs[bgcs/length(focal.cells)>0.001]
  bgcs <- bgcs[rev(order(bgcs))]
  pred <- BGC.pred
  values(X) <- NA
  values(X)[focal.cells] <- factor(pred, levels=units)[focal.cells]
  values(X)[1:length(units)] <- 1:length(units) # this is a patch that is necessary to get the color scheme right.
  
  plot(X, xaxt="n", yaxt="n", col=alpha(ColScheme, 1), legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE)
  plot(hill, col=grey(40:100/100), legend=FALSE, add=T)
  plot(X, add=T, xaxt="n", yaxt="n", col=alpha(ColScheme, 1), legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE)
  plot(focal.bdy, add=T, lwd=0.4)
  # plot(focal.mask, add=T, col=alpha("white", 1), border=F)
  mtext(GCM, side=3, line=-1.5, adj=0.99, cex=1, font=2)
  
  bgc.pct <- round(as.numeric(formatC(signif(bgcs/length(focal.cells)*100,digits=3), digits=3,format="fg", flag="#")),0)
  legendlength <- 10
  legendlength <- if(legendlength>length(bgcs)) length(bgcs) else legendlength  
  legend("bottomleft", cex=0.8, legend=paste(names(bgcs)[1:legendlength], " (", bgc.pct[1:legendlength], "%)", sep=""), fill=alpha(ColScheme[as.numeric(factor(names(bgcs), units))][1:legendlength], 1), bty="n")
  box()
  
  for(bgc in names(bgcs)){
    pts <- which(units[values(X)]==bgc)
    q=0.5
    pt <- xyFromCell(X, pts[min(which(pts >= quantile(pts, q)))])
    points(pt, pch=21, cex=1, lwd=0.8, bg=ColScheme[which(units==bgc)])
    text(pt-c(0, 0), bgc, pos=4, cex=0.6, font=2, offset=0.3)
    # print(q)
  }
  
  print(GCM)
}

dev.off()


################# hist year maps

png(filename=paste("Results\\CCISS", focal.unit, "BGCprojections.histyears.png",sep="."), type="cairo", units="in", width=9, height=6.5, pointsize=10, res=600)
par(mar=c(0.1,0.1, 0.1,0.1), mgp=c(2,0.25,0), mfrow=c(2,3))

for(hist.year in hist.years){

BGC.pred <- get(paste("BGC.pred", hist.year, sep="."))
zone.pred <- rep(NA, length(BGC.pred))
for(i in zones){ zone.pred[grep(i,BGC.pred)] <- i }
bgcs <- table(BGC.pred[focal.cells])
bgcs <- bgcs[bgcs/length(focal.cells)>0.001]
bgcs <- bgcs[rev(order(bgcs))]
pred <- BGC.pred
values(X) <- NA
values(X)[focal.cells] <- factor(pred, levels=units)[focal.cells]
values(X)[1:length(units)] <- 1:length(units) # this is a patch that is necessary to get the color scheme right.

plot(focal.bdy, lwd=0.4)
plot(hill, col=grey(40:100/100), legend=FALSE, add=T)
plot(X, add=T, xaxt="n", yaxt="n", col=alpha(ColScheme, 1), legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE)
plot(focal.bdy, add=T, lwd=0.4)
# plot(focal.mask, add=T, col=alpha("white", 1), border=F)
mtext(hist.year.name[which(hist.years==hist.year)], side=3, line=-1.5, adj=0.99, cex=1, font=2)

bgc.pct <- round(as.numeric(formatC(signif(bgcs/length(focal.cells)*100,digits=3), digits=3,format="fg", flag="#")),0)
legendlength <- 10
legendlength <- if(legendlength>length(bgcs)) length(bgcs) else legendlength  
legend("bottomleft", cex=1, legend=paste(names(bgcs)[1:legendlength], " (", bgc.pct[1:legendlength], "%)", sep=""), fill=alpha(ColScheme[as.numeric(factor(names(bgcs), units))][1:legendlength], 1), bty="n")
box()

for(bgc in names(bgcs)){
  pts <- which(units[values(X)]==bgc)
  q=0.5
  pt <- xyFromCell(X, pts[min(which(pts >= quantile(pts, q)))])
  points(pt, pch=21, cex=1, lwd=0.8, bg=ColScheme[which(units==bgc)])
  text(pt-c(0, 0), bgc, pos=4, cex=0.6, font=2, offset=0.3)
  # print(q)
}
print(hist.year)
}
dev.off()

#===============================================================================
# SUITABILITY -- variety plot for a single species
#===============================================================================
# spps.name <- c("Yellow cedar", "Interior spruce", "Ponderosa pine", "Lodgepole pine", "Western larch", "Western hemlock", "Mountain hemlock", "Douglas-fir", "Red alder", "Western redcedar", "Subalpine fir", "Grand fir")
spps.lookup <- read.csv("InputData\\Tree speciesand codes_2.0.csv")
edatope.name <- c("Subxeric-poor", "Mesic-medium", "Hygric-rich")
proj.year.name=c("2020s", "2050s", "2080s")
rcp=rcps[1]
proj.year=proj.years[2]

edatope="C4"
# for(edatope in edatopes){

spps <- c("Pl", "Fd", "Cw", "Ba", "Bl", "Bg", "Yc", "Pa", "Hm", "Lw", "Hw", "Py", "Dr", "Ep", "At", "Pw", "Ss", "Sb", "Qg", "Act")

# x11(width=6.5, height=8.5, pointsize=12)

spp=spps[1]
for(spp in spps){
  
  png(filename=paste("Results\\CCISS.manu.Suitability.basic",grid,spp, edatope, rcp, proj.year,"png",sep="."), type="cairo", units="in", width=6.5, height=6.7, pointsize=10, res=600)
  par(mar=c(0.1,0.1,0.1,0.1), mfrow=c(2,2), bg="white")
  
  RefSuit <- read.csv(paste("OutputData\\Suit.ref", grid, spp, edatope, "csv", sep="."))[,1]
  outRange.base <- RefSuit==5
  RefSuit[RefSuit==5] <- 4
  RefSuit[is.na(RefSuit)] <- 4
  
  # compile the GCM projections into a data frame
  ProjSuit <- data.frame(temp=rep(NA, length(RefSuit))) #initiate the data frame with a dummy column
  ChangeSuit <- data.frame(temp=rep(NA, length(RefSuit))) #initiate the data frame with a dummy column
  for(GCM in GCMs){
    temp <- read.csv(paste("OutputData\\Suit", grid, GCM, rcp, proj.year, spp, edatope, "csv", sep="."))
    temp[temp==5] <- 4
    temp[is.na(temp)] <- 4
    ProjSuit <- cbind(ProjSuit,temp)
    ChangeSuit <- cbind(ChangeSuit,RefSuit-temp)
  }
  ProjSuit <- ProjSuit[,-1] #remove the dummy column
  ChangeSuit <- ChangeSuit[,-1] #remove the dummy column
  names(ProjSuit) <- GCMs
  names(ChangeSuit) <- GCMs
  
  # calculate ensemble mean suitability. this isn't biased by missing suitabilties for exotic BGCs
  ChangeSuit.mean <- apply(ChangeSuit, 1, mean, na.rm=T)
  
  outRange <- outRange.base
  outRange[which(ChangeSuit.mean!=0)] <- FALSE
  ChangeSuit.mean[outRange==T] <- NA
  
  
  # map of historical suitability
  breakseq <- c(0.5,1.5,2.5,3.5,5)
  ColScheme <- c(brewer.pal(9,"Greys")[9], brewer.pal(9,"Greens")[7], brewer.pal(9,"Greens")[4], "white")
  ColScheme <- c("darkgreen", "dodgerblue1", "gold2", "white")
  length(ColScheme)
  
  values(X) <- NA
  values(X) <- RefSuit
  plot(X, xaxt="n", yaxt="n", col=alpha(ColScheme, 1), legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE)
  plot(focal.bdy, add=T, lwd=1.2)
  legend("topright", legend=c("1 (primary)", "2 (secondary)", "3 (tertiary)"), 
         fill=ColScheme, cex=0.9, title="Historical suitability", inset=0.015)
  # }
  # box()
  Common <- as.character(spps.lookup$EnglishName[which(spps.lookup$TreeCode==spp)])
  Latin <- as.character(spps.lookup$ScientificName[which(spps.lookup$TreeCode==spp)])
  panel <- paste("(", LETTERS[which(spps==spp)],")", sep="")
  mtext(if(spp%in%spps.lookup$TreeCode) bquote(bold(.(spp))~"-"~.(Common)) else bquote(.(panel)~bold(.(spp))),
        side=1, line=-1.75, adj=0.01, cex=0.8, font=2)
  # mtext(if(spp%in%spps.lookup$TreeCode) bquote(.(panel)~bold(.(spp))~"-"~.(Common)~"("*italic(.(Latin)*")")) else bquote(.(panel)~bold(.(spp))),
  #       side=3, line=-1.75, adj=0.01, cex=0.8, font=2)
  mtext(paste("Site type: ", edatope, " (", edatope.name[which(edatopes==edatope)], ")", sep=""), side=1, line=-2.75, adj=0.01, cex=0.7, font=1)
  box()
  
  # map of suitability change
  breakpoints <- seq(-3,3,0.5); length(breakpoints)
  labels <- c("-3","-2", "-1", "no change", "+1","+2","+3")
  ColScheme <- c(brewer.pal(11,"RdBu")[c(1,2,3,4,4)], "grey80", brewer.pal(11,"RdBu")[c(7,8,8,9,10,11)]); length(ColScheme)
  
  values(X)[focal.cells] <- ChangeSuit.mean[focal.cells]
  plot(focal.bdy, lwd=0.4)
  plot(X, add=T, xaxt="n", yaxt="n", col=alpha(ColScheme, 1), legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE)
  plot(focal.mask, add=T, border=F, col=alpha("white", 1))
  plot(focal.bdy, add=T, lwd=0.4)
  # if(spp==spps[1]){
  par(xpd=T)
  xl <- -121; yb <- 51.8; xr <- -120.85; yt <- 52.45
  rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme)
  text(rep(xr-.010000,length(labels)),seq(yb,yt,(yt-yb)/(length(labels)-1)),labels,pos=4,cex=0.8,font=1)
  text(xl-.050000, mean(c(yb,yt))-.030000, paste("Mean change\nin suitability (", proj.year.name[which(proj.years==proj.year)], ")", sep=""), srt=90, pos=3, cex=0.9, font=2)
  # }
  par(xpd=F)
  # box()
  
  
  # map of ensemble agreement on trend
  increasing <- which(ChangeSuit.mean>0)
  decreasing <- which(ChangeSuit.mean<0)
  bifurc <- rep(NA, length(ChangeSuit.mean))
  bifurc[outRange==F] <- 0
  bifurc[increasing] <- apply(ChangeSuit[increasing,], 1, function(x){return(sum(x< 0, na.rm=T)/length(x))})
  bifurc[decreasing] <- apply(ChangeSuit[decreasing,], 1, function(x){return(sum(x> 0, na.rm=T)/length(x))})
  values(X)[focal.cells] <- bifurc[focal.cells]
  
  breakpoints <- seq(0,0.5,0.1); length(breakpoints)
  labels <- c("High", "Medium", "Low")
  # ColScheme <- c(brewer.pal(11,"RdBu")[c(4,2,1)], "black", brewer.pal(11,"RdBu")[c(8,10,11)]); length(ColScheme)
  # ColScheme <- brewer.pal(11,"RdBu")[-c(2,5,7,10)]; length(ColScheme)
  ColScheme <- c("grey80", brewer.pal(9,"YlOrRd")[c(3,5,7,9)]); length(ColScheme)
  
  plot(focal.bdy, lwd=0.4)
  plot(X, add=T, xaxt="n", yaxt="n", col=alpha(ColScheme, 1), legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE)
  plot(focal.mask, add=T, border=F, col=alpha("white", 1))
  plot(focal.bdy, add=T, lwd=0.4)
  # mtext(paste("Edatope:", edatope), side=1, line=-1.5, adj=0.02, cex=1.1, font=2)
  # if(spp==spps[1]){
  xl <- -121; yb <- 51.8; xr <- -120.85; yt <- 52.45
  rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme)
  text(rep(xr,length(labels)),seq(yb,yt,(yt-yb)/(length(labels)*5-1))[c(3,13)],labels[c(1,3)],pos=4,cex=0.9,font=1)
  text(rep(xr,length(labels)),seq(yb,yt,(yt-yb)/(length(labels)*5-1))[c(1,15)],paste(round(breakpoints*100), "%", sep="")[c(1, length(breakpoints))],pos=4,cex=0.9,font=1)
  text(xl-.050000, mean(c(yb,yt))-.030000, paste("Ensemble agreement\non trend"), srt=90, pos=3, cex=1, font=2)
  # legend("bottomleft", legend=c(spp, paste("Edatope:", edatope), proj.year, rcp, " "), cex=1.4, bty="n", inset=-0.05)
  # }
  # box()
  
  # map of binary appearance/disappearance
  Suit.ensemble <- as.matrix(ProjSuit)
  Suit.ensemble[Suit.ensemble==5] <- 4
  binary <- rep(0, length(ProjSuit))
  binary[outRange.base==T] <- NA
  binary[outRange.base] <- apply(Suit.ensemble[outRange.base,], 1, function(x){return(if((sum(x<4, na.rm=T)/sum(!is.na(x)))>0) sum(x<4, na.rm=T)/sum(!is.na(x)) else NA)})
  binary[outRange.base==F] <- apply(Suit.ensemble[outRange.base==F,], 1, function(x){return(0-sum(x==4, na.rm=T)/sum(!is.na(x)))})
  values(X)[focal.cells] <- binary[focal.cells]
  
  breakpoints <- seq(-1,1,0.2); length(breakpoints)
  labels <- c("Retreat", "Expansion")
  ColScheme <- c(brewer.pal(11,"RdBu")[c(1:4)], "grey90", brewer.pal(11,"RdBu")[c(7:11)]); length(ColScheme)
  
  plot(focal.bdy, lwd=0.4)
  plot(X, add=T, xaxt="n", yaxt="n", col=alpha(ColScheme, 1), legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE)
  plot(focal.bdy, add=T, lwd=0.4)
  plot(focal.mask, add=T, border=F, col=alpha("white", 1))
  # mtext(paste("Edatope:", edatope), side=1, line=-1.5, adj=0.02, cex=1.1, font=2)
  # if(spp==spps[1]){
  xl <- -121; yb <- 51.8; xr <- -120.85; yt <- 52.45
  rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme)
  text(rep(xr+.010000,length(labels)),seq(yb,yt,(yt-yb)/(length(GCMs)-1))[c(2,9)],labels,pos=4,cex=0.9,font=0.8, srt=90)
  text(rep(xr-.020000,length(labels)),seq(yb,yt,(yt-yb)/(length(GCMs)-1))[c(1,8,15)],c("100%", "0%", "100%"),pos=4,cex=0.8,font=1)
  text(xl-.050000, mean(c(yb,yt))-.030000, paste("Change in presence/absence\n(% of models)"), srt=90, pos=3, cex=0.9, font=2)
  # legend("bottomleft", legend=c(spp, paste("Edatope:", edatope), proj.year, rcp, " "), cex=1.4, bty="n", inset=-0.05)
  # }
  # box()
  
  print(spp)
  dev.off()
  
}

par(mar=c(4,4,1,1), mfrow=c(1,1))
plot(dem)

#===============================================================================
# SUITABILITY -- Time series plot for a single species
#===============================================================================
# spps.name <- c("Yellow cedar", "Interior spruce", "Ponderosa pine", "Lodgepole pine", "Western larch", "Western hemlock", "Mountain hemlock", "Douglas-fir", "Red alder", "Western redcedar", "Subalpine fir", "Grand fir")
spps.lookup <- read.csv("InputData\\Tree speciesand codes_2.0.csv")
edatope.name <- c("Poor-Subxeric", "Medium-Mesic", "Rich-Hygric")
proj.year.name=c("2020s", "2050s", "2080s")
rcp=rcps[1]
proj.year=proj.years[2]

edatope="C4"
# for(edatope in edatopes){

spps <- c("Pl", "Fd", "Cw", "Ba", "Bl", "Bg", "Yc", "Pa", "Hm", "Lw", "Hw", "Py", "Dr", "Ep", "At", "Pw", "Ss", "Sb", "Qg", "Act")
spps <- c("Pl", "Fd", "Cw", "Bl", "Lw", "Hw", "Py", "Dr", "Ep", "At", "Pw", "Sb", "Act")

# x11(width=6.5, height=8.5, pointsize=12)

spp=spps[1]
for(edatope in edatopes){
  for(spp in spps){
    
    png(filename=paste("Results\\CCISS",focal.unit,"SuitabilityChange",spp, edatope,"png",sep="."), type="cairo", units="in", width=6.5, height=6.7, pointsize=10, res=600)
    par(mar=c(0.1,0.1,0.1,0.1), mfrow=c(2,2), bg="white")
    
    RefSuit <- read.csv(paste("OutputData\\Suit.ref", grid, spp, edatope, "csv", sep="."))[,1]
    outRange.base <- RefSuit==5
    RefSuit[RefSuit==5] <- 4
    RefSuit[is.na(RefSuit)] <- 4
    
    # map of historical suitability
    breakseq <- c(0.5,1.5,2.5,3.5,5)
    ColScheme <- c(brewer.pal(9,"Greys")[9], brewer.pal(9,"Greens")[7], brewer.pal(9,"Greens")[4], "white")
    ColScheme <- c("darkgreen", "dodgerblue1", "gold2", "white")
    length(ColScheme)
    
    values(X) <- NA
    values(X) <- RefSuit
    values(X)[1:4] <- 1:4 # this is a patch that is necessary to get the color scheme right.
    plot(X, xaxt="n", yaxt="n", col=alpha(ColScheme, 1), legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE)
    plot(focal.bdy, add=T, lwd=1.2)
    legend("topright", legend=c("1 (primary)", "2 (secondary)", "3 (tertiary)"), 
           fill=ColScheme, cex=0.9, title="Historical suitability", inset=0.015)
    # }
    box()
    
  
    # RECENT PERIOD
    hist.year <- 2009
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
    
    values(X) <- NA
    values(X) <- changeSuit
    plot(focal.bdy, lwd=0.4)
    plot(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks=breakpoints, legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE)
    plot(focal.mask, add=T, border=F, col=alpha("white", 1))
    plot(focal.bdy, add=T, lwd=0.4)
    # if(spp==spps[1]){
    par(xpd=T)
    xl <- -121; yb <- 51.8; xr <- -120.85; yt <- 52.45
    rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme)
    text(rep(xr,length(labels)),seq(yb,yt,(yt-yb)/(length(labels)*5-1))[c(2,8,14)],labels,pos=4,cex=1,font=1)
    text(xl-.050000, mean(c(yb,yt))-.030000, paste("Projected change\nin suitability (", hist.year.name[which(hist.years==hist.year)], ")", sep=""), srt=90, pos=3, cex=1, font=2)
    
    # }
    par(xpd=F)
    Common <- as.character(spps.lookup$EnglishName[which(spps.lookup$TreeCode==spp)])
    Latin <- as.character(spps.lookup$ScientificName[which(spps.lookup$TreeCode==spp)])
    panel <- paste("(", LETTERS[which(spps==spp)],")", sep="")
    mtext(if(spp%in%spps.lookup$TreeCode) bquote(bold(.(spp))~"-"~.(Common)) else bquote(.(panel)~bold(.(spp))),
          side=1, line=-2.5, adj=0.01, cex=1, font=2)
    # mtext(if(spp%in%spps.lookup$TreeCode) bquote(.(panel)~bold(.(spp))~"-"~.(Common)~"("*italic(.(Latin)*")")) else bquote(.(panel)~bold(.(spp))),
    #       side=3, line=-1.75, adj=0.01, cex=0.8, font=2)
    mtext(paste("Site type: ", edatope, " (", edatope.name[which(edatopes==edatope)], ")", sep=""), side=1, line=-1.25, adj=0.01, cex=0.9, font=1)
    
    
    # FUTURE PERIODS
    for(proj.year in c(2025, 2055)){
      # compile the GCM projections into a data frame
      ProjSuit <- data.frame(temp=rep(NA, length(RefSuit))) #initiate the data frame with a dummy column
      ChangeSuit <- data.frame(temp=rep(NA, length(RefSuit))) #initiate the data frame with a dummy column
      for(GCM in GCMs){
        temp <- read.csv(paste("OutputData\\Suit", grid, GCM, rcp, proj.year, spp, edatope, "csv", sep="."))[,1]
        temp[temp==5] <- 4
        temp[is.na(temp)] <- 4
        ProjSuit <- cbind(ProjSuit,temp)
        ChangeSuit <- cbind(ChangeSuit,RefSuit-temp)
      }
      ProjSuit <- ProjSuit[,-1] #remove the dummy column
      ChangeSuit <- ChangeSuit[,-1] #remove the dummy column
      names(ProjSuit) <- GCMs
      names(ChangeSuit) <- GCMs
      
      # calculate ensemble mean suitability. this isn't biased by missing suitabilties for exotic BGCs
      ChangeSuit.mean <- apply(ChangeSuit, 1, mean, na.rm=T)
      
      outRange <- outRange.base
      outRange[which(ChangeSuit.mean!=0)] <- FALSE
      ChangeSuit.mean[outRange==T] <- NA
      
      # map of suitability change
      breakpoints <- seq(-3,3,0.5); length(breakpoints)
      labels <- c("-3","-2", "-1", "no change", "+1","+2","+3")
      ColScheme <- c(brewer.pal(11,"RdBu")[c(1,2,3,4,4)], "grey80", brewer.pal(11,"RdBu")[c(7,8,8,9,10,11)]); length(ColScheme)
      
      values(X) <- ChangeSuit.mean
      plot(focal.bdy, lwd=0.4)
      plot(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks=breakpoints, legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE)
      plot(focal.mask, add=T, border=F, col=alpha("white", 1))
      plot(focal.bdy, add=T, lwd=0.4)
      # if(spp==spps[1]){
      par(xpd=T)
      xl <- -121; yb <- 51.8; xr <- -120.85; yt <- 52.45
      rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme)
      text(rep(xr-.010000,length(labels)),seq(yb,yt,(yt-yb)/(length(labels)-1)),labels,pos=4,cex=0.8,font=1)
      text(xl-.050000, mean(c(yb,yt))-.030000, paste("Mean change\nin suitability (", proj.year.name[which(proj.years==proj.year)], ")", sep=""), srt=90, pos=3, cex=0.9, font=2)
      # }
      par(xpd=F)
      # box()
    }
    
    dev.off()
    print(spp)
  }
  print(edatope)
}


