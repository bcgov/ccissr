##======================================================================================
## CCISS Publication Scripts
## Step 2 - Species suitability projections
##======================================================================================

# Colin Mahony
# c_mahony@alumni.ubc.ca
# 778-288-4008
# July 21, 2019

# What this script does: calculate suitability outputs (but not stocking standards) equivalent to the CCISS tool, but simplified for fast running. 
# reason for script: the full CCISS script takes ~10-20 seconds per POI, meaning weeks of running time for a publishable raster grid. 
# Objective: produce results for a ~250,000-cell grid at a run time of maximum 8 hours.  
# general approach: this script achieves speed by: 
# 1. running only for selected edatopes
# 2. running on the whole grid vector, not looping on the POIs. 
# 3. outputs and inputs in many small tables, rather than a few big tables. 
# 4. running only for selected species

# Analysis Steps: 
#   /Edatope/RCP/Spp/proj.year/GCM: 
#     1. Project BGC unit 
#     2. find the site series associated with the edatope
#     3. obtain suitability of spp in site series

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


# rm(list=ls())

#===============================================================================
# Set analysis Parameters
#===============================================================================

##setwd("C:\\Colin\\Projects\\2019_CCISS")
setwd("C:/Users/kirid/Desktop/2019_CCISS")

grid <- "BC2kmGrid"

GCMs <-  c("ACCESS1-0","CanESM2","CCSM4","CESM1-CAM5","CNRM-CM5","CSIRO-Mk3-6-0", "GFDL-CM3","GISS-E2R", "HadGEM2-ES", "INM-CM4", "IPSL-CM5A-MR", "MIROC-ESM", "MIROC5", "MPI-ESM-LR","MRI-CGCM3")
rcps <- c("rcp45", "rcp85")
proj.years <- c(2025, 2055, 2085)
hist.years <- c(1995, 2004, 2005, 2009, 2014, 2017)
edatopes<- c("B2", "C4", "D6")
edatope.name <- c("Subxeric-poor", "Mesic-medium", "Hygric-rich")
BGCcolors <- read.csv("C:\\Colin\\Projects\\2019_CCISS\\InputData\\BGCzone_Colorscheme.csv")

# Knowledge Tables
inputFolder = "inputs"
treesuit="TreeSpp_ESuit_v11_18"
SiteSeries_Use <-read.csv(paste(inputFolder,"SiteSeries_Use_5.csv",sep="/"),stringsAsFactors=FALSE,na.strings=".")
spps.lookup <- read.csv(paste(inputFolder,"Tree speciesand codes_2.0_2May2019.csv", sep = "/"))


#===============================================================================
# generate the vector of mapped BGCs
#===============================================================================

points <- fread(paste(inputFolder,"/", grid,".csv", sep=""))
BGC <- points$ID2
BGC <- gsub(" ","",BGC)  
sort(table(BGC))

#===============================================================================
# Import BGC projections for each period
#===============================================================================

## reference period BGC
BGC.pred.ref <- as.character(read.csv(paste("outputs\\BGC.pred.ref", grid, "csv", sep="."))[,1])

# Historical BGC
for(hist.year in hist.years){
  BGC.pred <- as.character(read.csv(paste("outputs\\BGC.pred", grid,hist.year,"csv", sep="."))[,1])
  assign(paste("BGC.pred",hist.year,sep="."), BGC.pred)
  print(hist.year)
}

# Future BGC
predDat <- foreach(rcp = rcps, .combine = rbind) %:%
  foreach(proj.year = proj.years, .combine = rbind) %:%
  foreach(GCM = GCMs, .combine = rbind) %do% {
    BGC.pred <- fread(paste("outputs\\BGC.pred",grid, GCM, rcp, proj.year,".csv", sep=""))
    BGC.pred <- BGC.pred %>% mutate(Scn = rcp, FuturePeriod = proj.year, GCM = GCM)
    cat(GCM, rcp, proj.year, "\n")
    BGC.pred
  }


#===============================================================================
# import will's lookup table for the site series associated with the selected edatope in each BGC unit
#===============================================================================

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

hist.years <- c(1995, 2004, 2005, 2009, 2014, 2017)

# Import suitability tables
S1 <- read.csv(paste("inputs/",treesuit,".csv",sep=""),stringsAsFactors=F,na.strings=".")
S1 <- unique(S1)

## EDA: are there suitabilities for all projected units? 
# NoSuit <- PredSum[-which(PredSum$BGC%in%S1$BGC),]
# NoSuit[rev(order(NoSuit$count)),]

## EDA: Which site series are missing suitabilities? 
for(edatope in edatopes) assign(paste("NoSuit", edatope, sep="."), SiteLookup[-which(SiteLookup[,which(names(SiteLookup)==edatope)]%in%S1$Unit),which(names(SiteLookup)==edatope)])
for(edatope in edatopes) print(get(paste("NoSuit", edatope, sep=".")))

## EDA: are there any units missing from the SiteSeries_Use table? 
BGClist <- unique(S1$BGC)
BGClist[-which(BGClist%in%SiteLookup$MergedBGC)]

# select the species to run the analysis on
spps <- unique(S1$Spp)
spps <- spps[-which(spps=="X")]
spps.candidate <- spps.lookup$TreeCode[-which(spps.lookup$Exclude=="x")]
spps <- spps[which(spps%in%spps.candidate)] 

require(doParallel)
library(iterators)
set.seed(123321)
coreNum <- as.numeric(detectCores()-2)
cl <- makeCluster(coreNum)
registerDoParallel(cl, cores = coreNum)

outAll <- foreach(spp = spps,.combine = rbind) %do% {
  foreach(edatope = edatopes, .combine = rbind) %do% {
    BGC.pred <- as.character(BGC.pred.ref) # get the BGC prediction
    
    # get the suitability for the selected species associated with each site series
    suit <- S1$ESuit[which(S1$Spp==spp)][match(as.vector(unlist(SiteLookup[which(names(SiteLookup)==edatope)])), S1$Unit[which(S1$Spp==spp)])]
    Suit.ref <- suit[match(BGC.pred, SiteLookup$MergedBGC)]
    Suit.ref[is.na(Suit.ref)] <- 5 #set the NA values to suitability 5
    Suit.ref[Suit.ref==4] <- 5 #set 4 to suitability 5
    ##fwrite(Suit.ref, paste("outputs\\Suit.ref", grid, spp, edatope, "csv", sep="."), row.names = F)
    
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
      write.csv(temp, paste("outputs\\Suit", grid, hist.year, spp, edatope, "csv", sep="."), row.names = F)
      # print(hist.year)
    }
    
    cat("starting",spp,edatope, "\n")

    # get the suitability for future periods, for each projection/edatope/species combination
    suit <- foreach(dat = isplit(predDat$BGC.pred, list(predDat$GCM,predDat$Scn,predDat$FuturePeriod)), .combine = rbind,.packages = c("randomForest","dplyr", "foreach","iterators")) %dopar% {
          # get the BGC projection and sub in the crosswalk between the modeled units and the table units
          BGC.pred <- dat$value
          # BGC.pred[which(BGC.pred%in%Crosswalk$Modeled)] <- as.character(Crosswalk$Tables[match(BGC.pred[which(BGC.pred%in%Crosswalk$Modeled)], Crosswalk$Modeled)]) # XXX THIS IS NOT CORRECT. NEED TO FIGURE OUT HOW TO INCORPORATE THE CROSSWALK TABLE PROPERLY. sub in the crosswalk between the modeled units and the table units
          
          ## identify cells with no suitability interpretations
          bgc.exotic <- (1:length(BGC.pred))[-which(BGC.pred%in%unique(BGC))]
          bgc.exotic.noSuit <- bgc.exotic[-which(BGC.pred[bgc.exotic] %in% unique(S1$BGC))]
          
          # get the suitability for the selected species associated with each site series
          suit <- S1$ESuit[which(S1$Spp==spp)][match(as.vector(unlist(SiteLookup[which(names(SiteLookup)==edatope)])), S1$Unit[which(S1$Spp==spp)])]
          temp <- suit[match(BGC.pred, SiteLookup$MergedBGC)]
          temp[is.na(temp)] <- 5 #set the NA values to suitability 5 (weights unsuitable a bit more heavily than suitable classes during averaging)
          temp[temp==4] <- 5 #set 4 to suitability 5
          temp[bgc.exotic.noSuit] <- NA # set cells with no suitabilty interpretatoin to NA
          out <- data.frame(Suit = temp) %>% mutate(GCM = dat[["key"]][[1]], Scn = dat[["key"]][[2]], FuturePeriod = dat[["key"]][[3]], Edatope = edatope, Spp = spp)
          ##cat(unlist(dat$key), "\n")
          out
    }
  }
}

fwrite(outAll, "BGC.pred_AllFuture.csv")
