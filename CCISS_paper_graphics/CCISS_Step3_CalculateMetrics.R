
##======================================================================================
## CCISS Publication Scripts
## Step 3 - calculation/export of community metrics (suitability richness, turnover and persistence) and scatterplots against MAT change. 
##======================================================================================

# Colin Mahony
# c_mahony@alumni.ubc.ca
# 778-288-4008
# July 21, 2019


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
library(msir)

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

# Knowledge Tables
treesuit="TreeSpp_ESuit_v11_18"
SiteSeries_Use <-read.csv(paste("InputData/","SiteSeries_Use_5",".csv",sep=""),stringsAsFactors=FALSE,na.strings=".")
spps.lookup <- read.csv("InputData\\Tree speciesand codes_2.0_2May2019.csv")

#===============================================================================
# generate the vector of mapped BGCs
#===============================================================================

points <- read.csv(paste("InputData\\",grid,".csv", sep=""))
BGC <- points$ID2
table(BGC)

BGC <- gsub(" ","",BGC)  
table(BGC)
sort(table(BGC))

# # merge small units
# BGC[which(BGC=="CMAwh")] <- "CMAun"
# BGC[which(BGC=="MHun")] <- "MHunp"
# BGC[which(BGC=="ESSFxvw")] <- "ESSFxvp"
# BGC[which(BGC=="ESSFdcp")] <- "ESSFdcw"

#BGC zones
BGCcolors <- read.csv("C:\\Colin\\Projects\\2019_CCISS\\InputData\\BGCzone_Colorscheme.csv")
zone <- rep(NA, length(BGC))
for(i in BGCcolors$zone){ zone[grep(i,BGC)] <- i }
table(zone)

#===============================================================================
# Import suitability tables
#===============================================================================
wd="InputData"
treesuit2=paste(wd,"/",treesuit,".csv",sep="")
S1 <- read.csv(treesuit2,stringsAsFactors=F,na.strings=".")
S1 <- unique(S1)

# select the species to run the analysis on
spps <- unique(S1$Spp)
spps <- spps[-which(spps=="X")]
spps.candidate <- spps.lookup$TreeCode[-which(spps.lookup$Exclude=="x")]
spps <- spps[which(spps%in%spps.candidate)] 


#===============================================================================
# assemble the species mix for the reference period
#===============================================================================
for(edatope in edatopes){
  comm.ref <- as.data.frame(matrix(rep(NA, dim(points)[1]*length(spps)), dim(points)[1], length(spps)))
  for(spp in spps){
    Suit <- read.csv(paste("OutputData\\Suit.ref", grid, spp, edatope, "csv", sep="."))[,1]
    Suit[is.na(Suit)] <- 5  #XXX note this is different from the equivalent line for the other time periods. 
    Suit <- 1-(Suit-1)/4
    comm.ref[,which(spps==spp)] <- Suit
  }
  names(comm.ref) <- spps
  comm.ref.spp <- comm.ref>0
  
  assign(paste("comm.ref", edatope, sep=""), comm.ref)

  SuitRichness <- apply(comm.ref, 1, sum, na.rm=T)
  SppRichness <- apply(comm.ref.spp, 1, sum, na.rm=T)
  write.csv(SuitRichness, paste("OutputData\\SuitRichness.ref", grid, edatope, "csv", sep="."), row.names = F)
  write.csv(SppRichness, paste("OutputData\\SppRichness.ref", grid, edatope, "csv", sep="."), row.names = F)

  print(edatope)
}

#===============================================================================
# for each hist.year, assemble the community and calculate richness and turnover. write to disk
#===============================================================================
# hist.years <- c(2004, 2009)

for(edatope in edatopes){
  comm.ref <- get(paste("comm.ref", edatope, sep=""))
  comm.ref.spp <- comm.ref>0

  for(hist.year in hist.years){

    comm.hist <- as.data.frame(matrix(rep(NA, dim(points)[1]*length(spps)), dim(points)[1], length(spps)))
    for(spp in spps){
      Suit <- read.csv(paste("OutputData\\Suit", grid, hist.year, spp, edatope, "csv", sep="."))[,1]
      Suit[Suit==4] <- 5
      Suit <- 1-(Suit-1)/4
      comm.hist[,which(spps==spp)] <- Suit
    }
    names(comm.hist) <- spps
    comm.hist.spp <- comm.hist>0
    
    SuitRichness <- apply(comm.hist, 1, sum, na.rm=T)
    SppRichness <- apply(comm.hist.spp, 1, sum, na.rm=T)
    write.csv(SuitRichness, paste("OutputData\\SuitRichness", grid, hist.year, edatope, "csv", sep="."), row.names = F)
    write.csv(SppRichness, paste("OutputData\\SppRichness", grid, hist.year, edatope, "csv", sep="."), row.names = F)

    #suitability turnover
    SuitTurnover <- apply(abs(comm.hist-comm.ref), 1, sum, na.rm=T)/apply(cbind(comm.ref,comm.hist), 1, sum, na.rm=T)
    SuitTurnover[which(apply(comm.ref, 1, sum, na.rm=T)==0)] <- 99
    SuitTurnover[which(apply(cbind(comm.ref,comm.hist), 1, sum, na.rm=T)==0)] <- NA
    write.csv(SuitTurnover, paste("OutputData\\SuitTurnover", grid, hist.year, edatope, "csv", sep="."), row.names = F)

    #species turnover (presence/absence)
    SppTurnover <- apply(abs(comm.hist.spp-comm.ref.spp), 1, sum, na.rm=T)/apply(cbind(comm.ref.spp,comm.hist.spp), 1, sum, na.rm=T)
    SppTurnover[which(apply(comm.ref.spp, 1, sum, na.rm=T)==0)] <- 99
    SppTurnover[which(apply(cbind(comm.ref.spp,comm.hist.spp), 1, sum, na.rm=T)==0)] <- NA
    write.csv(SppTurnover, paste("OutputData\\SppTurnover", grid, hist.year, edatope, "csv", sep="."), row.names = F)

    #suitability persistence
    comm <- cbind(comm.hist, comm.ref)
    SuitPersistence <- apply(comm, 1, FUN=function(x){return(sum(x[which(x[(length(x)/2+1):length(x)]>0)])/sum(x[(length(x)/2+1):length(x)]))})
    SuitPersistence[which(apply(comm.ref, 1, sum, na.rm=T)==0)] <- NA
    write.csv(SuitPersistence, paste("OutputData\\SuitPersistence", grid, hist.year, edatope, "csv", sep="."), row.names = F)
    
    #species persistence
    comm <- cbind(comm.hist.spp, comm.ref.spp)
    SppPersistence <- apply(comm, 1, FUN=function(x){return(sum(x[which(x[(length(x)/2+1):length(x)]>0)])/sum(x[(length(x)/2+1):length(x)]))})
    SppPersistence[which(apply(comm.ref.spp, 1, sum, na.rm=T)==0)] <- NA
    write.csv(SppPersistence, paste("OutputData\\SppPersistence", grid, hist.year, edatope, "csv", sep="."), row.names = F)
    
    print(hist.year)
  }
  print(edatope)
}

# ===============================================================================
# for each GCM/edatope/rcp/proj.year, assemble the community and calculate richness and turnover. write to disk
# ===============================================================================


edatope <- "C4"
rcp <- "rcp45"
proj.year <- 2055


for(GCM in GCMs){
  for(edatope in edatopes){
    comm.ref <- get(paste("comm.ref", edatope, sep=""))
    comm.ref.spp <- comm.ref>0
    
    for(rcp in rcps){
      for(proj.year in proj.years){

        comm.proj <- as.data.frame(matrix(rep(NA, dim(points)[1]*length(spps)), dim(points)[1], length(spps)))
        for(spp in spps){
          Suit <- read.csv(paste("OutputData\\Suit", grid, GCM, rcp, proj.year, spp, edatope, "csv", sep="."))[,1]
          Suit[Suit==4] <- 5
          Suit <- 1-(Suit-1)/4
          comm.proj[,which(spps==spp)] <- Suit
        }
        names(comm.proj) <- spps
        comm.proj.spp <- comm.proj>0
        
        SuitRichness <- apply(comm.proj, 1, sum, na.rm=T)
        SppRichness <- apply(comm.proj.spp, 1, sum, na.rm=T)
        write.csv(SuitRichness, paste("OutputData\\SuitRichness", grid, GCM, rcp, proj.year, edatope, "csv", sep="."), row.names = F)
        write.csv(SppRichness, paste("OutputData\\SppRichness", grid, GCM, rcp, proj.year, edatope, "csv", sep="."), row.names = F)
        
        #suitability turnover
        SuitTurnover <- apply(abs(comm.proj-comm.ref), 1, sum, na.rm=T)/apply(cbind(comm.ref,comm.proj), 1, sum, na.rm=T)
        SuitTurnover[which(apply(comm.ref, 1, sum, na.rm=T)==0)] <- 99
        SuitTurnover[which(apply(cbind(comm.ref,comm.proj), 1, sum, na.rm=T)==0)] <- NA
        write.csv(SuitTurnover, paste("OutputData\\SuitTurnover", grid, GCM, rcp, proj.year, edatope, "csv", sep="."), row.names = F)

        #species turnover (presence/absence)
         SppTurnover <- apply(abs(comm.proj.spp-comm.ref.spp), 1, sum, na.rm=T)/apply(cbind(comm.ref.spp,comm.proj.spp), 1, sum, na.rm=T)
        SppTurnover[which(apply(comm.ref.spp, 1, sum, na.rm=T)==0)] <- 99
        SppTurnover[which(apply(cbind(comm.ref.spp,comm.proj.spp), 1, sum, na.rm=T)==0)] <- NA
        write.csv(SppTurnover, paste("OutputData\\SppTurnover", grid, GCM, rcp, proj.year, edatope, "csv", sep="."), row.names = F)
        
        #suitability persistence
        comm <- cbind(comm.proj, comm.ref)
        SuitPersistence <- apply(comm, 1, FUN=function(x){return(sum(x[which(x[(length(x)/2+1):length(x)]>0)])/sum(x[(length(x)/2+1):length(x)]))})
        SuitPersistence[which(apply(comm.ref, 1, sum, na.rm=T)==0)] <- NA
        write.csv(SuitPersistence, paste("OutputData\\SuitPersistence", grid, GCM, rcp, proj.year, edatope, "csv", sep="."), row.names = F)
        
        #species persistence
        comm <- cbind(comm.proj.spp, comm.ref.spp)
        SppPersistence <- apply(comm, 1, FUN=function(x){return(sum(x[which(x[(length(x)/2+1):length(x)]>0)])/sum(x[(length(x)/2+1):length(x)]))})
        SppPersistence[which(apply(comm.ref.spp, 1, sum, na.rm=T)==0)] <- NA
        write.csv(SppPersistence, paste("OutputData\\SppPersistence", grid, GCM, rcp, proj.year, edatope, "csv", sep="."), row.names = F)
        
        # print(proj.year)
      }
      # print(rcp)
    }
    print(edatope)
  }
  print(GCM)
}

## Note that species turnover can either be lower or higher than suitability turnover.
par(mar=c(4,4,1,1))
sample <- sample(1:length(BGC), 1000)
plot(SuitTurnover[sample], SppTurnover[sample], xlim=c(0,1), ylim=c(0,1))
lines(c(0,1), c(0,1))

## Note that there is no relationship between richness and turnover
richness <- get(paste("richness.ref", edatope, sep=""))
plot(SuitTurnover[sample], richness[sample], xlim=c(0,1))
plot(SppTurnover[sample], richness[sample], xlim=c(0,1))

## Note the strong relationship between richness and richness change
richness.change <- apply(comm.proj, 1, sum)-richness
plot(richness[sample], richness.change[sample])
## This is improved somewhat by expressing change proportionally
richness.changePct <- apply(comm.proj, 1, sum)/richness
plot(richness[sample], log(richness.changePct[sample]))
lines(c(-99,99), c(0,0))


#===============================================================================
# calculate mean MAT change for each model prediction
#===============================================================================

fplot=paste("InputData\\", grid, "_Normal_1961_1990MSY.csv", sep="")
Y0 <- fread(fplot, select = "MAT", stringsAsFactors = FALSE, data.table = FALSE) #fread is faster than read.csv
MAT.ref <- Y0$MAT
MAT.mean.ref <- mean(MAT.ref, na.rm=T)

for(hist.year in hist.years){
  Y0 <- fread(paste("InputData\\", grid, "_", hist.year, "_BioVars.csv", sep=""), select = "MAT", stringsAsFactors = FALSE, data.table = FALSE) #fread is faster than read.csv
  assign(paste("MAT", hist.year, sep="."), Y0$MAT)
 assign(paste("MAT.change", hist.year, sep="."), mean(Y0$MAT, na.rm=T)-MAT.mean.ref)
  print(hist.year)
}

for(GCM in GCMs){
  Y1 <- fread(paste("InputData\\", grid, "_", GCM, "_BioVars.csv", sep=""), select = c("GCM", "MAT"), stringsAsFactors = FALSE, data.table = FALSE)
  ## assign single vectors to RCPs and proj.years
  Ystr <- strsplit(Y1[,1], "_")
  Y4 <- matrix(unlist(Ystr), ncol=3, byrow=TRUE)
  for(rcp in rcps){
    for(proj.year in proj.years){
      assign(paste("MAT", GCM, rcp, proj.year, sep="."), Y1$MAT[which(Y4[,2]==rcp & Y4[,3]==proj.year)])
    }
    print(rcp)
  }
  print(GCM)
}

# calculate mean climate values for each GCM/year/rcp
for(rcp in rcps){
  for(proj.year in proj.years){
    MAT.change.mean <- rep(NA, length(GCMs))
    for(GCM in GCMs){
      MAT <- get(paste("MAT", GCM, rcp, proj.year, sep="."))
      MAT.change <- MAT-MAT.mean.ref
      MAT.change.mean[which(GCMs==GCM)] <- mean(MAT.change, na.rm=T)
    }
    assign(paste("MAT.change", rcp, proj.year, sep="."), MAT.change.mean)
    print(proj.year)
  }
  print(rcp)
}

#===============================================================================
# calculate mean suitability turnover for each GCM/year/rcp
#===============================================================================

for(edatope in edatopes){
  for(hist.year in hist.years){
    SuitTurnover <- read.csv(paste("OutputData\\SuitTurnover", grid, hist.year, edatope, "csv", sep="."))[,1]
    SuitTurnover[SuitTurnover==99] <- NA # remove grid cells with no current suitability
    SuitTurnover.mean <- mean(SuitTurnover, na.rm=T)
    SppTurnover <- read.csv(paste("OutputData\\SppTurnover", grid, hist.year, edatope, "csv", sep="."))[,1]
    SppTurnover[SppTurnover==99] <- NA # remove grid cells with no current suitability
    SppTurnover.mean <- mean(SppTurnover, na.rm=T)
    assign(paste("SuitTurnover.mean", hist.year, edatope, sep="."), SuitTurnover.mean)
    assign(paste("SppTurnover.mean", hist.year, edatope, sep="."), SppTurnover.mean)
  }
  print(edatope)
}


for(edatope in edatopes){
  for(rcp in rcps){
    for(proj.year in proj.years){
      SuitTurnover.mean <- rep(NA, length(GCMs))
      SppTurnover.mean <- rep(NA, length(GCMs))
      for(GCM in GCMs){
        SuitTurnover <- read.csv(paste("OutputData\\SuitTurnover", grid, GCM, rcp, proj.year, edatope, "csv", sep="."))[,1]
        SuitTurnover[SuitTurnover==99] <- NA # remove grid cells with no current suitability
                SuitTurnover.mean[which(GCMs==GCM)] <- mean(SuitTurnover, na.rm=T)
        SppTurnover <- read.csv(paste("OutputData\\SppTurnover", grid, GCM, rcp, proj.year, edatope, "csv", sep="."))[,1]
        SppTurnover[SppTurnover==99] <- NA # remove grid cells with no current suitability
        SppTurnover.mean[which(GCMs==GCM)] <- mean(SppTurnover, na.rm=T)
      }
      assign(paste("SuitTurnover.mean", rcp, proj.year, edatope, sep="."), SuitTurnover.mean)
      assign(paste("SppTurnover.mean", rcp, proj.year, edatope, sep="."), SppTurnover.mean)
      print(proj.year)
    }
    print(rcp)
  }
  print(edatope)
}



#===============================================================================
# Suitability turnover relative to temperature change (predicted baseline)
#===============================================================================
hist.periods <- c("1991-2000", "1991-2017", "2001-2010", "2001-2017","2011-2017", "2017")

edatope=edatopes[2]
for(edatope in edatopes){
png(filename=paste("Results\\CCISS.SuitTurnoverVsMATchange",edatope, "png",sep="."), type="cairo", units="in", width=3.5, height=3.5, pointsize=9, res=300)
par(mfrow=c(1,1), mar=c(3.25,3.25,0.1,0.1), mgp=c(2.25,0.25,0))
plot(0, xlim=c(0,7.2), ylim=c(0,.85), yaxs="i", xaxs="i", col="white", xaxt="n", yaxt="n", 
     xlab=bquote(BC~mean~temperature~change~relative~to~"1961-90"~"("*degree*C*")"), 
     ylab="BC mean suitability turnover")
axis(1, at=0:8, labels = 0:8, tck=0)
axis(2, at=seq(0,1,0.2), labels = paste(seq(0,1,0.2)*100, "%", sep=""), las=2, tck=0)

for(rcp in rcps){
  for(proj.year in proj.years){
    temp.suit <- rep(NA, length(GCMs))
    temp.spp <- rep(NA, length(GCMs))
    SuitTurnover.mean <- get(paste("SuitTurnover.mean", rcp, proj.year, edatope, sep=".")) 
    SppTurnover.mean <- get(paste("SppTurnover.mean", rcp, proj.year, edatope, sep=".")) 
    
    x <- get(paste("MAT.change", rcp, proj.year, sep="."))
    y1 <- SuitTurnover.mean
    # y1 <- SppTurnover.mean
    
    points(x,y1, pch=c(21,22)[which(rcps==rcp)], bg=c("black", "dodgerblue", "yellow")[which(proj.years==proj.year)])
    # points(x,y2, col="red", pch=16)
    # points(0,mean(BGC.pred.ref!=BGC), pch=16, cex=1.3)
    # points(0,mean(zone.pred.ref!=zone), col="red", pch=16, cex=1.3)
  }
}

for(hist.year in hist.years[c(1,2,3,4,5,6)]){
  x <- get(paste("MAT.change", hist.year, sep="."))
  y <- get(paste("SuitTurnover.mean", hist.year, edatope, sep="."))
  # y <- get(paste("SppTurnover.mean", hist.year, edatope, sep="."))
  points(x,y, cex=1.4, pch=2, col="grey")
  text(x,y, hist.periods[which(hist.years==hist.year)], pos=4, cex=0.8, font=2, col="grey")
  print(y)
  }

legend("bottomright", legend=c("RCP4.5", "RCP8.5", "2011-2040", "2041-2070", "2071-2100"), pch=c(21,22, NA,NA,NA), 
       pt.bg=c("gray", "gray", NA,NA,NA), pt.cex=1.5, fill=c(NA, NA, "black", "dodgerblue", "yellow"), border = F, bty="n")

dev.off()
}



#===============================================================================
# Species turnover relative to temperature change (predicted baseline)
#===============================================================================
hist.periods <- c("1991-2000", "1991-2017", "2001-2010", "2001-2017","2011-2017", "2017")

edatope=edatopes[2]
for(edatope in edatopes){
  png(filename=paste("Results\\CCISS.SppTurnoverVsMATchange",edatope, "png",sep="."), type="cairo", units="in", width=3.5, height=3.5, pointsize=9, res=300)
  par(mfrow=c(1,1), mar=c(3.25,3.25,0.1,0.1), mgp=c(2.25,0.25,0))
  plot(0, xlim=c(0,7.2), ylim=c(0,.85), yaxs="i", xaxs="i", col="white", xaxt="n", yaxt="n", 
       xlab=bquote(BC~mean~temperature~change~relative~to~"1961-90"~"("*degree*C*")"), 
       ylab="BC mean species turnover")
  axis(1, at=0:8, labels = 0:8, tck=0)
  axis(2, at=seq(0,1,0.2), labels = paste(seq(0,1,0.2)*100, "%", sep=""), las=2, tck=0)
  
  for(rcp in rcps){
    for(proj.year in proj.years){
      temp.suit <- rep(NA, length(GCMs))
      temp.spp <- rep(NA, length(GCMs))
      SuitTurnover.mean <- get(paste("SuitTurnover.mean", rcp, proj.year, edatope, sep=".")) 
      SppTurnover.mean <- get(paste("SppTurnover.mean", rcp, proj.year, edatope, sep=".")) 
      
      x <- get(paste("MAT.change", rcp, proj.year, sep="."))
      # y1 <- SuitTurnover.mean
      y1 <- SppTurnover.mean
      
      points(x,y1, pch=c(21,22)[which(rcps==rcp)], bg=c("black", "dodgerblue", "yellow")[which(proj.years==proj.year)])
      # points(x,y2, col="red", pch=16)
      # points(0,mean(BGC.pred.ref!=BGC), pch=16, cex=1.3)
      # points(0,mean(zone.pred.ref!=zone), col="red", pch=16, cex=1.3)
    }
  }
  
  for(hist.year in hist.years[c(1,2,3,4,5,6)]){
    x <- get(paste("MAT.change", hist.year, sep="."))
    # y <- get(paste("SuitTurnover.mean", hist.year, edatope, sep="."))
    y <- get(paste("SppTurnover.mean", hist.year, edatope, sep="."))
    points(x,y, cex=1.4, pch=2, col="grey")
    text(x,y, hist.periods[which(hist.years==hist.year)], pos=4, cex=0.8, font=2, col="grey")
    print(y)
  }
  
  legend("bottomright", legend=c("RCP4.5", "RCP8.5", "2011-2040", "2041-2070", "2071-2100"), pch=c(21,22, NA,NA,NA), 
         pt.bg=c("gray", "gray", NA,NA,NA), pt.cex=1.5, fill=c(NA, NA, "black", "dodgerblue", "yellow"), border = F, bty="n")
  
  dev.off()
}



#===============================================================================
# calculate mean persistence for each GCM/year/rcp
#===============================================================================

for(edatope in edatopes){
  for(hist.year in hist.years){
    SuitPersistence <- read.csv(paste("OutputData\\SuitPersistence", grid, hist.year, edatope, "csv", sep="."))[,1]
    SuitPersistence[SuitPersistence==99] <- NA # remove grid cells with no current suitability
    SuitPersistence.mean <- mean(SuitPersistence, na.rm=T)
    SppPersistence <- read.csv(paste("OutputData\\SppPersistence", grid, hist.year, edatope, "csv", sep="."))[,1]
    SppPersistence[SppPersistence==99] <- NA # remove grid cells with no current suitability
    SppPersistence.mean <- mean(SppPersistence, na.rm=T)
    assign(paste("SuitPersistence.mean", hist.year, edatope, sep="."), SuitPersistence.mean)
    assign(paste("SppPersistence.mean", hist.year, edatope, sep="."), SppPersistence.mean)
  }
  print(edatope)
}


for(edatope in edatopes){
  for(rcp in rcps){
    for(proj.year in proj.years){
      SuitPersistence.mean <- rep(NA, length(GCMs))
      SppPersistence.mean <- rep(NA, length(GCMs))
      for(GCM in GCMs){
        SuitPersistence <- read.csv(paste("OutputData\\SuitPersistence", grid, GCM, rcp, proj.year, edatope, "csv", sep="."))[,1]
        SuitPersistence[SuitPersistence==99] <- NA # remove grid cells with no current suitability
        SuitPersistence.mean[which(GCMs==GCM)] <- mean(SuitPersistence, na.rm=T)
        SppPersistence <- read.csv(paste("OutputData\\SppPersistence", grid, GCM, rcp, proj.year, edatope, "csv", sep="."))[,1]
        SppPersistence[SppPersistence==99] <- NA # remove grid cells with no current suitability
        SppPersistence.mean[which(GCMs==GCM)] <- mean(SppPersistence, na.rm=T)
      }
      assign(paste("SuitPersistence.mean", rcp, proj.year, edatope, sep="."), SuitPersistence.mean)
      assign(paste("SppPersistence.mean", rcp, proj.year, edatope, sep="."), SppPersistence.mean)
      print(proj.year)
    }
    print(rcp)
  }
  print(edatope)
}

#vectors of projection specs. 
seq.rcp <- NA
seq.proj.year <-  NA
seq.GCM <-  NA
for(rcp in rcps){
  for(proj.year in proj.years){
    for(GCM in GCMs){
      seq.rcp <- c(seq.rcp, rcp)
      seq.proj.year <- c(seq.proj.year, proj.year)
      seq.GCM <- c(seq.GCM, GCM)
    }
  }
}

#===============================================================================
# Suitability Persistence relative to temperature change (predicted baseline)
#===============================================================================
hist.periods <- c("1991-2000", "1991-2017", "2001-2010", "2001-2017","2011-2017", "2017")

edatope=edatopes[2]
for(edatope in edatopes){
  png(filename=paste("Results\\CCISS.SuitPersistenceVsMATchange",edatope, "png",sep="."), type="cairo", units="in", width=3.5, height=3.5, pointsize=9, res=300)
  par(mfrow=c(1,1), mar=c(3.25,3.25,0.1,0.1), mgp=c(2.25,0.25,0))
  plot(0, xlim=c(0,7.2), ylim=c(0,1.2), yaxs="i", xaxs="i", col="white", xaxt="n", yaxt="n", 
       xlab=bquote(BC~mean~temperature~change~relative~to~"1961-90"~"("*degree*C*")"), 
       ylab="BC mean suitability Persistence")
  axis(1, at=0:8, labels = 0:8, tck=0)
  axis(2, at=seq(0,1.2,0.2), labels = paste(seq(0,1.2,0.2)*100, "%", sep=""), las=2, tck=0)
  
  for(rcp in rcps){
    for(proj.year in proj.years){
      temp.suit <- rep(NA, length(GCMs))
      temp.spp <- rep(NA, length(GCMs))
      SuitPersistence.mean <- get(paste("SuitPersistence.mean", rcp, proj.year, edatope, sep=".")) 
      SppPersistence.mean <- get(paste("SppPersistence.mean", rcp, proj.year, edatope, sep=".")) 
      
      x <- get(paste("MAT.change", rcp, proj.year, sep="."))
      y1 <- SuitPersistence.mean
      # y1 <- SppPersistence.mean
      
      points(x,y1, pch=c(21,22)[which(rcps==rcp)], bg=c("black", "dodgerblue", "yellow")[which(proj.years==proj.year)])
      # points(x,y2, col="red", pch=16)
      # points(0,mean(BGC.pred.ref!=BGC), pch=16, cex=1.3)
      # points(0,mean(zone.pred.ref!=zone), col="red", pch=16, cex=1.3)
    }
  }
  
  for(hist.year in hist.years[c(1,2,3,4,5,6)]){
    x <- get(paste("MAT.change", hist.year, sep="."))
    y <- get(paste("SuitPersistence.mean", hist.year, edatope, sep="."))
    # y <- get(paste("SppPersistence.mean", hist.year, edatope, sep="."))
    points(x,y, cex=1.4, pch=2, col="grey")
    text(x,y, hist.periods[which(hist.years==hist.year)], pos=4, cex=0.8, font=2, col="grey")
    print(y)
  }
  
  legend("bottomright", legend=c("RCP4.5", "RCP8.5", "2011-2040", "2041-2070", "2071-2100"), pch=c(21,22, NA,NA,NA), 
         pt.bg=c("gray", "gray", NA,NA,NA), pt.cex=1.5, fill=c(NA, NA, "black", "dodgerblue", "yellow"), border = F, bty="n")
  
  dev.off()
}

#===============================================================================
# Species Persistence relative to temperature change (predicted baseline)
#===============================================================================
hist.periods <- c("1991-2000", "1991-2017", "2001-2010", "2001-2017","2011-2017", "2017")

edatope=edatopes[2]
for(edatope in edatopes){
  png(filename=paste("Results\\CCISS.SppPersistenceVsMATchange",edatope, "png",sep="."), type="cairo", units="in", width=3.5, height=3.5, pointsize=9, res=300)
  par(mfrow=c(1,1), mar=c(3.25,3.25,0.1,0.1), mgp=c(2.25,0.25,0))
  plot(0, xlim=c(0,7.2), ylim=c(0,1.2), yaxs="i", xaxs="i", col="white", xaxt="n", yaxt="n", 
       xlab=bquote(BC~mean~temperature~change~relative~to~"1961-90"~"("*degree*C*")"), 
       ylab="BC mean species persistence")
  axis(1, at=0:8, labels = 0:8, tck=0)
  axis(2, at=seq(0,1.2,0.2), labels = paste(seq(0,1.2,0.2)*100, "%", sep=""), las=2, tck=0)
  
  for(rcp in rcps){
    for(proj.year in proj.years){
      temp.suit <- rep(NA, length(GCMs))
      temp.spp <- rep(NA, length(GCMs))
      SuitPersistence.mean <- get(paste("SuitPersistence.mean", rcp, proj.year, edatope, sep=".")) 
      SppPersistence.mean <- get(paste("SppPersistence.mean", rcp, proj.year, edatope, sep=".")) 
      
      x <- get(paste("MAT.change", rcp, proj.year, sep="."))
      # y1 <- SuitPersistence.mean
      y1 <- SppPersistence.mean
      
      points(x,y1, pch=c(21,22)[which(rcps==rcp)], bg=c("black", "dodgerblue", "yellow")[which(proj.years==proj.year)])
      # points(x,y2, col="red", pch=16)
      # points(0,mean(BGC.pred.ref!=BGC), pch=16, cex=1.3)
      # points(0,mean(zone.pred.ref!=zone), col="red", pch=16, cex=1.3)
    }
  }
  
  for(hist.year in hist.years[c(1,2,3,4,5,6)]){
    x <- get(paste("MAT.change", hist.year, sep="."))
    # y <- get(paste("SuitPersistence.mean", hist.year, edatope, sep="."))
    y <- get(paste("SppPersistence.mean", hist.year, edatope, sep="."))
    points(x,y, cex=1.4, pch=2, col="grey")
    text(x,y, hist.periods[which(hist.years==hist.year)], pos=4, cex=0.8, font=2, col="grey")
    print(y)
  }
  
  legend("bottomright", legend=c("RCP4.5", "RCP8.5", "2011-2040", "2041-2070", "2071-2100"), pch=c(21,22, NA,NA,NA), 
         pt.bg=c("gray", "gray", NA,NA,NA), pt.cex=1.5, fill=c(NA, NA, "black", "dodgerblue", "yellow"), border = F, bty="n")
  
  dev.off()
}


#===============================================================================
# Species Persistence relative to temperature change (predicted baseline)
#===============================================================================
hist.periods <- c("1991-2000", "1991-2017", "2001-2010", "2001-2017","2011-2017", "2017")

metrics <- c("SuitPersistence", "SppPersistence")
metric=metrics[1]
edatope=edatopes[2]
# for(metric in metrics){
  # png(filename=paste("Results\\CCISS.MATchangeVs",metric, "png",sep="."), type="cairo", units="in", width=3.5, height=3.5, pointsize=9, res=300)
  par(mfrow=c(1,1), mar=c(3.25,3.25,0.1,0.1), mgp=c(2.25,0.25,0))
  plot(0, xlim=c(0,7.2), ylim=c(0,1.2), yaxs="i", xaxs="i", col="white", xaxt="n", yaxt="n", 
       xlab=bquote(BC~mean~temperature~change~relative~to~"1961-90"~"("*degree*C*")"), 
       ylab="BC mean species Persistence")
  axis(1, at=0:8, labels = 0:8, tck=0)
  axis(2, at=seq(0,1.2,0.2), labels = paste(seq(0,1.2,0.2)*100, "%", sep=""), las=2, tck=0)
  
  for(edatope in edatopes){
    temp.suit <- rep(NA, length(GCMs))
    temp.spp <- rep(NA, length(GCMs))
    for(rcp in rcps){
      for(proj.year in proj.years){
        temp.suit <- rep(NA, length(GCMs))
        temp.spp <- rep(NA, length(GCMs))
        
        x <- get(paste("MAT.change", rcp, proj.year, sep="."))
        y <- get(paste(metric, "mean", rcp, proj.year, edatope, sep=".")) 
        
        l <- loess.sd(y~x, span=1, nsigma = 1)
        polygon(c(l$x, rev(l$x)), c(l$upper, rev(l$lower)), col=alpha((1:3)[which(edatopes==edatope)], 0.5), border=NA)
        l <- loess(y~x)
        par(xpd=T)
        text(max(x)-0.25, predict(l, max(x)), edatope, pos=4, font= if(edatope=="C4") 2 else 1, col=ColScheme[which(edatopes==edatope)], cex= if(edatope=="C4") 1.1 else 1)
        par(xpd=F)

        if(edatope=="C4") points(x,y, pch=c(21,22)[which(rcps==rcp)], bg=c("black", "dodgerblue", "yellow")[which(proj.years==proj.year)])
        # points(x,y2, col="red", pch=16)
        # points(0,mean(BGC.pred.ref!=BGC), pch=16, cex=1.3)
        # points(0,mean(zone.pred.ref!=zone), col="red", pch=16, cex=1.3)
        
        
              }
    }
  }
  
  for(hist.year in hist.years[c(2,4)]){
    x <- get(paste("MAT.change", hist.year, sep="."))
    # y <- get(paste("SuitPersistence.mean", hist.year, edatope, sep="."))
    y <- get(paste(metric, "mean", hist.year, edatope, sep="."))
    points(x,y, cex=1.4, pch=2)
    text(x,y, hist.periods[which(hist.years==hist.year)], pos=4, cex=0.8, font=2)
    print(y)
  }
  
  legend("bottomright", legend=c("RCP4.5", "RCP8.5", "2011-2040", "2041-2070", "2071-2100"), pch=c(21,22, NA,NA,NA), 
         pt.bg=c("gray", "gray", NA,NA,NA), pt.cex=1.5, fill=c(NA, NA, "black", "dodgerblue", "yellow"), border = F, bty="n")
  
  # dev.off()
# }





