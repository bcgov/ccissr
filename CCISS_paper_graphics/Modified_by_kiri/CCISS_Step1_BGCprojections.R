
## ======================================================================================
## CCISS Publication Scripts Step 1 - BGC projections for BC study area
## ======================================================================================

# Colin Mahony c_mahony@alumni.ubc.ca 778-288-4008 July 21, 2019

require(RGtk2)
require(plyr)
require(rChoiceDialogs)
require(data.table)
require(doBy)
require(utils)
require(labdsv)
require(tools)
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
require(data.table)


rm(list = ls())

# ===============================================================================
# Set analysis Parameters
# ===============================================================================

grid <- "BC2kmGrid"

GCMs <- c("ACCESS1-0", "CanESM2", "CCSM4", "CESM1-CAM5", "CNRM-CM5", "CSIRO-Mk3-6-0", 
          "GFDL-CM3", "GISS-E2R", "HadGEM2-ES", "INM-CM4", "IPSL-CM5A-MR", "MIROC-ESM", 
          "MIROC5", "MPI-ESM-LR", "MRI-CGCM3")

addVars <- function(dat){
  dat$PPT_MJ <- dat$PPT05 + dat$PPT06  # MaY/June precip
  dat$PPT_JAS <- dat$PPT07 + dat$PPT08 + dat$PPT09  # July/Aug/Sept precip
  dat$PPT.dormant <- dat$PPT_at + dat$PPT_wt  # for calculating spring deficit
  dat$CMD.def <- 500 - (dat$PPT.dormant)  # start of growing season deficit original value was 400 but 500 seems better
  dat$CMD.def[dat$CMD.def < 0] <- 0  #negative values set to zero = no deficit
  dat$CMDMax <- dat$CMD07
  dat$CMD.total <- dat$CMD.def + dat$CMD
  return(dat)
}


### Load random forest model
model <- "5.1"
fname <- "inputs/BGCv11_AB_USA_16VAR_SubZone_RFmodel.Rdata"
load(fname)
#rownames(importance(BGCmodel)) ### shows the variable used in the RFmodel


GCMs <- c("ACCESS1-0", "CanESM2", "CCSM4", "CESM1-CAM5", "CNRM-CM5", "CSIRO-Mk3-6-0", 
          "GFDL-CM3", "GISS-E2R", "HadGEM2-ES", "INM-CM4", "IPSL-CM5A-MR", "MIROC-ESM", 
          "MIROC5", "MPI-ESM-LR", "MRI-CGCM3")
rcps <- c("rcp45", "rcp85")
proj.years <- c(2025, 2055, 2085)
hist.years <- c(1995, 2004, 2005, 2009, 2014, 2017)
edatopes <- c("B2", "C4", "D6")
spps.lookup <- fread("inputs/Tree speciesand codes_2.0_2May2019.csv")
edatope.name <- c("Subxeric-poor", "Mesic-medium", "Hygric-rich")
BGCcolors <- fread("inputs/BGCzone_Colorscheme.csv")


# ===============================================================================
# generate the vector of mapped BGCs
# ===============================================================================

points <- fread(paste("inputs/", grid, ".csv", sep = ""))
BGC <- points$ID2
table(BGC)

BGC <- gsub(" ", "", BGC)
table(BGC)
sort(table(BGC))

# # merge small units BGC[which(BGC=='CMAwh')] <- 'CMAun'
# BGC[which(BGC=='MHun')] <- 'MHunp' BGC[which(BGC=='ESSFxvw')] <-
# 'ESSFxvp' BGC[which(BGC=='ESSFdcp')] <- 'ESSFdcw'

# BGC zones

zone <- rep(NA, length(BGC))
for (i in BGCcolors$zone){
  zone[grep(i, BGC)] <- i
}
table(zone)


# #===============================================================================
# # select variables from ClimateWNA file and subset by GCM # NB: is
# commented out because you only need to do this once: it just reduces
# ram needed to load in the ClimateBC data
# #===============================================================================
# Columns <- unique(c('PPT05', 'PPT06', 'PPT07', 'PPT08', 'PPT09', 'PPT_at',
# 'PPT_wt', 'CMD07', 'CMD', 'MAT', 'PPT_sm', 'Tmin_wt', 'Tmax_sm',
# rownames(importance(BGCmodel))[-which(rownames(importance(BGCmodel))%in%c('PPT_MJ',
# 'PPT_JAS', 'PPT.dormant', 'CMD.def', 'CMDMax', 'CMD.total'))]))
# 
# #first batch of 8 models fplot=paste('inputs/', grid,
# '_48GCMs_MSYT.csv', sep='') Y1 <- fread(fplot, select = c('GCM',
# Columns), stringsAsFactors = FALSE, data.table = FALSE) #fread is
# faster than read.csv models <-
# c('ACCESS1-0','CanESM2','CCSM4','CESM1-CAM5','CNRM-CM5','CSIRO-Mk3-6-0',
# 'GFDL-CM3','GISS-E2R') for(model in models){ temp <-
# Y1[grep(model,Y1$GCM),] write.csv(temp, paste('inputs/', grid, '_',
# model, '_BioVars.csv', sep=''), row.names = F)
# print(which(models==model)) } #second batch of 7 models
# fplot=paste('inputs/', grid, '_42GCMs_MSYT.csv', sep='') Y1 <-
# fread(fplot, select = c('GCM', Columns), stringsAsFactors = FALSE,
# data.table = FALSE) #fread is faster than read.csv models <-
# c('HadGEM2-ES', 'INM-CM4', 'IPSL-CM5A-MR', 'MIROC-ESM', 'MIROC5',
# 'MPI-ESM-LR','MRI-CGCM3') for(model in models){ temp <-
# Y1[grep(model,Y1$GCM),] write.csv(temp, paste('inputs/', grid, '_',
# model, '_BioVars.csv', sep=''), row.names = F)
# print(which(models==model)) }


# ===============================================================================
# BGC Projections for reference period
# ===============================================================================

# setwd('C:/GitHub/2019_CCISS')
Columns <- unique(c("PPT05", "PPT06", "PPT07", "PPT08", "PPT09", "PPT_at", 
                    "PPT_wt", "CMD07", "CMD", "MAT", "PPT_sm", "Tmin_wt", "Tmax_sm",
                    rownames(importance(BGCmodel))[-which(rownames(importance(BGCmodel)) %in%
                                                            c("PPT_MJ", "PPT_JAS", "PPT.dormant", "CMD.def", "CMDMax", "CMD.total"))]))

fplot <- paste("inputs/", grid, "_Normal_1961_1990MSY.csv", sep = "")

Y0 <- fread(fplot, select = Columns, stringsAsFactors = FALSE, data.table = FALSE)  #fread is faster than read.csv

Y0 <- Y0[!is.na(Y0[, 2]), ]

Y0 <- addVars(Y0)

## Predict future subzones######
BGC.pred.ref <- predict(BGCmodel, Y0)
dir.create("./outputs")
write.csv(BGC.pred.ref, paste("outputs/BGC.pred.ref", grid, "csv", sep = "."), 
          row.names = F)

## Write Climate file ######
fwrite(Y0, paste("inputs/", grid, "_1961_1990_BioVars.csv", sep = ""))

# ===============================================================================
# BGC Projections for historical decades
# ===============================================================================

# setwd('C:\\Colin\\Projects\\2019_CCISS')
hist.years <- c(1995, 2005)
hist.periods <- c("1991_2000", "2001_2010")

for (hist.year in hist.years){
  hist.period <- hist.periods[which(hist.years == hist.year)]
  fplot <- paste("inputs/", grid, "_Decade_", hist.period, "MSY.csv", 
                 sep = "")
  
  Y0 <- fread(fplot, select = Columns, stringsAsFactors = FALSE, data.table = FALSE)  #fread is faster than read.csv
  Y0 <- Y0[!is.na(Y0[, 2]), ]
  # str(Y0)
  Y0 <- addVars(Y0)
  
  ## Predict future subzones######
  assign(paste("BGC.pred", hist.year, sep = "."), predict(BGCmodel, Y0))
  write.csv(get(paste("BGC.pred", hist.year, sep = ".")), paste("outputs/BGC.pred", 
                                                                grid, hist.year, "csv", sep = "."), row.names = F)
  
  ## Write Climate file ######
  write.csv(Y0, paste("inputs/", grid, "_", hist.year, "_BioVars.csv", 
                      sep = ""), row.names = F)
  
  print(hist.year)
}

# ===============================================================================
# BGC Projections for last decade
# ===============================================================================

# setwd('C:\\Colin\\Projects\\2019_CCISS')
fplot <- paste("inputs/", grid, "_2011-2018MSY.csv", sep = "")

Y0 <- fread(fplot, select = c("ID1", "Year", Columns), stringsAsFactors = FALSE, 
            data.table = FALSE)  #fread is faster than read.csv
# Y0 <- Y0[!is.na(Y0[,2]),]
str(Y0)

Y0 <- addVars(Y0)

# Extract the year 2017
Y2017 <- Y0[which(Y0$Year == 2017), ]
str(Y2017)

# Calculate mean of 2011-2017 period
Y1 <- Y0 %>%
  group_by(ID1) %>%
  summarise_all(list(mean)) %>%
  ungroup()
##Y1 <- aggregate(Y0, by = list(Y0$ID1), FUN = mean, na.rm = T)[, -1]
Y1 <- Y1[match(Y2017$ID1, Y1$ID1), ]
Y1 <- Y1[,-2]

## Predict BGC units######
BGC.pred.2017 <- predict(BGCmodel, Y2017)
BGC.pred.2014 <- predict(BGCmodel, Y1)
write.csv(BGC.pred.2017, paste("outputs/BGC.pred", grid, "2017.csv", sep = "."), 
          row.names = F)
write.csv(BGC.pred.2014, paste("outputs/BGC.pred", grid, "2014.csv", sep = "."), 
          row.names = F)

## Write Climate file ######
fwrite(Y2017, paste("inputs/", grid, "_2017_BioVars.csv", sep = ""))
fwrite(Y1, paste("inputs/", grid, "_2014_BioVars.csv", sep = ""))

# ===============================================================================
# BGC Projections for other historical normals
# ===============================================================================

# setwd('C:\\Colin\\Projects\\2019_CCISS')
hist.years <- c(1995, 2005)
hist.periods <- c("1991_2000", "2001_2010")

mean(c(1991, 2017))
mean(c(2001, 2017))

# read in the data for the 1990s and 2000s
for (hist.year in hist.years){
  hist.period <- hist.periods[which(hist.years == hist.year)]
  fplot <- paste("inputs/", grid, "_Decade_", hist.period, "MSY.csv", 
                 sep = "")
  
  Y0 <- fread(fplot, select = Columns, stringsAsFactors = FALSE, data.table = FALSE)  #fread is faster than read.csv
  Y0 <- Y0[!is.na(Y0[, 2]), ]
  
  Y0 <- addVars(Y0)
  assign(paste("Y", hist.year, sep = "."), Y0)
  print(hist.year)
}

# 2011-2017 period already exported in previous phase, so read that in
Y.2014 <- read.csv(paste("inputs/", grid, "_2014_BioVars.csv", sep = ""))
Y.2014 <- Y.2014[, -which(names(Y.2014) %in% c("ID1", "Year"))]
str(Y.2014)

Y.2009 <- (Y.2005 * 10 + Y.2014 * 7)/17
Y.2004 <- (Y.1995 * 10 + Y.2005 * 10 + Y.2014 * 7)/27
str(Y.2009)
str(Y.2004)

hist.years <- c(2004, 2009)
for (hist.year in hist.years){
  
  ## Predict future subzones######
  assign(paste("BGC.pred", hist.year, sep = "."), predict(BGCmodel, get(paste("Y", 
                                                                              hist.year, sep = "."))))
  write.csv(get(paste("BGC.pred", hist.year, sep = ".")), paste("outputs/BGC.pred", 
                                                                grid, hist.year, "csv", sep = "."), row.names = F)
  
  ## Write Climate file ######
  write.csv(get(paste("Y", hist.year, sep = ".")), paste("inputs/", grid, 
                                                         "_", hist.year, "_BioVars.csv", sep = ""), row.names = F)
  
  print(hist.year)
}

# ===============================================================================
# BGC Projections for future periods
# ===============================================================================

library(tidyr)
Y0 <- fread(paste("inputs/",grid,"_90 GCMsMSY.csv", sep = ""),select = c("Year", "ID1","ID2", Columns))   ##
Y0 <- separate(Y0, Year, into = c("Model","Scenario","FuturePeriod"), sep = "_", remove = T)
Y0$FuturePeriod <- gsub(".gcm","",Y0$FuturePeriod)

require(doParallel)
set.seed(123321)
coreNum <- as.numeric(detectCores()-2)
cl <- makeCluster(coreNum)
registerDoParallel(cl, cores = coreNum)

out <- foreach(GCM = GCMs, .combine = rbind) %:%
  foreach(rcp = rcps, .combine = rbind) %:%
  foreach(proj.year = proj.years, .combine = rbind, .packages = c("data.table","randomForest", "dplyr")) %do% {
    sub <- Y0[Y0$Model == GCM & Y0$Scenario == rcp & Y0$FuturePeriod == proj.year,]
    sub <- addVars(sub)
    subPred <- data.frame(ID = sub$ID1)
    subPred$BGC.pred <- predict(BGCmodel, sub)
    fwrite(subPred, paste("outputs/BGC.pred", grid, GCM, rcp, proj.year, ".csv", sep = ""))
    temp <- aggregate(ID ~ BGC.pred, data = subPred, FUN = length) %>% 
      mutate(Model = GCM, Scn = rcp, FuturePeriod = proj.year)
    temp
  }



