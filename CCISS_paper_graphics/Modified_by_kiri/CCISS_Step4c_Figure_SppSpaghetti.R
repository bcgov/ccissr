
##======================================================================================
## CCISS Publication Scripts
## Step 4c - Figure - Spaghetti plots of spp suitable area relative to MAT change
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
edatope.name <- c("nutrient poor, subxeric", "nutrient-medium, mesic", "nutrient-rich, hygric")
proj.year.name=c("2020s", "2050s", "2080s")
rcp.name=c("RCP4.5", "RCP8.5")

## non-THLB BGCs for exclusion from results
points <- read.csv(paste("InputData\\",grid,".csv", sep=""))
BGC <- points$ID2
BGC <- gsub(" ","",BGC)
BGCs_notin_THLB <- read.csv("InputData\\BGCs_notin_THLB.csv")
BGCs_notin_THLB <- BGCs_notin_THLB$BGC[which(BGCs_notin_THLB$Exlude=="x")]
exclude <- which(BGC%in%BGCs_notin_THLB)

spps.lookup <- read.csv("InputData\\Tree speciesand codes_2.0_2May2019.csv")

#===============================================================================
# Import suitability tables
#===============================================================================
wd="InputData"
treesuit="TreeSpp_ESuit_v11_18"
treesuit2=paste(wd,"/",treesuit,".csv",sep="")
S1 <- read.csv(treesuit2,stringsAsFactors=F,na.strings=".")
S1 <- unique(S1)

# select a subset of species to run the analysis on
spps <- unique(S1$Spp)
spps.candidate <- spps.lookup$TreeCode[-which(spps.lookup$Exclude=="x")]
spps <- spps[which(spps%in%spps.candidate)]



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

# Compile the MAT for all time periods/scenarios into a single vector
MAT.change <- vector()
for(rcp in rcps){
  for(proj.year in proj.years){
    MAT.change <- c(MAT.change, get(paste("MAT.change", rcp, proj.year, sep=".")))
    # print(proj.year)
  }
  # print(rcp)
}

############################
#calculate the total number of grid cells occupied by each species in each projection

for(spp in spps){
  for(edatope in edatopes){
    Suit <- read.csv(paste("OutputData\\Suit.ref", grid, spp, edatope, "csv", sep="."))[-exclude,1]
    Suit[is.na(Suit)] <- 5
    Suit <- 1-(Suit-1)/4
    assign(paste("SuitCells.ref", spp, edatope, sep="."), sum(Suit))
    assign(paste("SppCells.ref", spp, edatope, sep="."), sum(Suit>0))
    # print(edatope)
  }
  print(spp)
}


for(spp in spps){
  for(edatope in edatopes){
    for(rcp in rcps){
      for(proj.year in proj.years){
        
        SuitCells <- rep(NA, length(GCMs))
        SppCells <- rep(NA, length(GCMs))
        for(GCM in GCMs){
          
          Suit.proj <- read.csv(paste("OutputData\\Suit", grid, GCM, rcp, proj.year, spp, edatope, "csv", sep="."))[-exclude,1]
          Suit.proj[is.na(Suit.proj)] <- 5
          Suit.proj <- 1-(Suit.proj-1)/4
          
          SuitCells[which(GCMs==GCM)] <- sum(Suit.proj)
          SppCells[which(GCMs==GCM)] <- sum(Suit.proj>0)
          # print(GCM)
        }
        
        assign(paste("SuitCells", rcp, proj.year, spp, edatope, sep="."), SuitCells)
        assign(paste("SppCells", rcp, proj.year, spp, edatope, sep="."), SppCells)
        # print(proj.year)
      }
      # print(rcp)
    }
    # print(edatope)
  }
  print(paste(spp, " ",round(which(spps==spp)/length(spps)*100,1), "%", sep=""))
}

# Compile the metrics for all time periods/scenarios into a single vector
for(edatope in edatopes){
  for(spp in spps){
    SuitCells <- vector()
    SppCells <- vector()
    for(rcp in rcps){
      for(proj.year in proj.years){
        SuitCells <- c(SuitCells, get(paste("SuitCells", rcp, proj.year, spp, edatope, sep=".")))
        SppCells <- c(SppCells, get(paste("SppCells", rcp, proj.year, spp, edatope, sep=".")))
        # print(proj.year)
      }
      # print(rcp)
    }
    assign(paste("SuitCells", spp, edatope, sep="."), SuitCells)
    assign(paste("SppCells", spp, edatope, sep="."), SppCells)
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


############################
## analysis of exotic species suitable area
############################
rcp.focal <- "rcp4.5"
proj.year.focal <- 2085

length(BGC[-exclude])

exotic.table <- data.frame(spp=NA, 
                           edatope=NA,
                           area.pct=NA, 
                           EnglishName=NA, 
                           ScientificName=NA
                           )
for(edatope in edatopes){
  exotic.area.pct <- vector()
spps.exotic <- spps[-which(spps%in%spps.native)]
for(spp in spps.exotic){
  temp <- get(paste("SuitCells", spp, edatope, sep="."))[which(seq.rcp==rcp.focal & seq.proj.year==proj.year.focal)]
  exotic.area.pct[which(spps.exotic==spp)] <- mean(temp)/length(BGC[-exclude])
}
temp <- data.frame(spp=spps.exotic[rev(order(exotic.area.pct))], 
                   edatope=rep(edatope, length(spps.exotic)),
           area.pct=paste(round(exotic.area.pct[rev(order(exotic.area.pct))]*100,2), "%", sep=""), 
           EnglishName=spps.lookup$EnglishName[match(spps.exotic[rev(order(exotic.area.pct))], spps.lookup$TreeCode)], 
           ScientificName=spps.lookup$ScientificName[match(spps.exotic[rev(order(exotic.area.pct))], spps.lookup$TreeCode)])
exotic.table <- rbind(exotic.table, temp)
print(edatope)
}
exotic.table <- exotic.table[rev(order(exotic.table$area.pct)),]




############################
## two-panel plot of species trends relative to MAT change, log and raw scaled
############################
native <- spps.lookup$Native[match(spps,spps.lookup$TreeCode)]
spps.native <- c("Ra", "Pl", "Pj", "Fd", "Cw", "Ba", "Sx", "Bl", "Bg", "Yc", "Pa", "Hm", "Lw", "La", "Lt", "Hw", "Py", "Dr", "Ep", "At", "Acb", "Pw", "Ss", "Sb", "Qg", "Act", "Mb")

for(edatope in edatopes){

  png(filename=paste("Results\\CCISS_manu_SppSpaghetti", edatope, "png",sep="."), type="cairo", units="in", width=6.5, height=5, pointsize=8, res=400)
  
  for(spp in spps){
    x <- c(0,MAT.change)
    y <- c(get(paste("SuitCells.ref", spp, edatope, sep=".")),get(paste("SuitCells", spp, edatope, sep=".")))*4 #times 4km^2 because they are 2km grid cells. 
    l <- loess(y[order(x)]~x[order(x)])
    assign(paste("line",spp, sep="."), predict(l, seq(0,max(MAT.change), 0.01)))
  }
  
  par(mfrow=c(1,2))
  for(transform in c(T, F)){
    
    par(mar=c(5,4,0,0), mgp=c(4, 0.2, 0))
    ylim=if(transform==T) c(2,6.1) else c(0,47000)
    plot(0, xlim=c(-2.5,9.5), ylim=ylim, yaxs="i", xaxs="i", col="white", xaxt="n", yaxt="n", 
         xlab=bquote(BC~mean~temperature~change~relative~to~"1961-90"~"("*degree*C*")"), 
         ylab="")
    axis(1, at=0:7, labels = 0:7, tck=0)
    par(mgp=c(3, 0.2, 0))
    if(transform==T) title(ylab="Suitable area (sq. km)")
    axis(2, at=if(transform==T) 0:7 else seq(0,50000,10000), labels = if(transform==T) format(10^(0:7), scientific = FALSE, big.mark=",") else format(seq(0,50000,10000), scientific = FALSE, big.mark=","), tck=0, las=2)
    # rect(-9,0,0, 60000, col="lightgray", border=F)
    # rect(max(MAT.change),0,9, ylim[2]*1.1, col="lightgray", border=F)
    
    suit.exotic.final <- vector()
    for(spp in spps[-which(spps%in%spps.native)]){
      line <- get(paste("line",spp, sep="."))
      if(transform==T) line[line<1] <- 1
      if(transform==T) line <- log10(line)
      suit.exotic.final[which(spps[-which(spps%in%spps.native)]==spp)] <- line[length(line)]
    }
    
    suit.native.initial <- vector()
    for(spp in spps[which(spps%in%spps.native)]){
      suit.native.initial[which(spps[which(spps%in%spps.native)]==spp)] <- get(paste("SuitCells.ref", spp, edatope, sep="."))
    }
    
    spplist <- spps[which(spps%in%spps.native)][order(suit.native.initial)]
    colors = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)][-1]
    set.seed(2)
    ColScheme <- c(brewer.pal(n=12, "Paired"),sample(colors,length(spps)-12))
    for(spp in spplist){
      i <- which(spplist==spp)
      line <- get(paste("line",spp, sep="."))
      if(transform==T) line[line<1] <- 1
      if(transform==T) line <- log10(line)
      if(line[1]> if(transform==T) ylim[1] else 100){
        lines(seq(0,max(MAT.change), 0.01), line, col=ColScheme[i], lwd=2)
        position <- rep(0:3, times=100)
        text(0-position[i]*0.6, line[1], spp, pos=2, col=ColScheme[i], font=2, cex=0.75, offset=0.1)
        lines(c(0-position[i]*0.6,0), rep(line[1],2), col=ColScheme[i], lty=2)
      }
    }
    
    spplist <- spps[-which(spps%in%spps.native)][order(suit.exotic.final)]
    for(spp in spplist){
      i <- which(spplist==spp)
      line <- get(paste("line",spp, sep="."))
      if(transform==T) line[line<1] <- 1
      if(transform==T) line <- log10(line)
      if(max(line)> if(transform==T) ylim[1] else 100){
        lines(seq(0,max(MAT.change), 0.01), line)
        position <- rep(0:3, times=100)
        if(which.max(line)>(length(line)-50)){
          text(max(MAT.change)+position[i]*0.6, line[length(line)], spp, pos=4, cex=0.75, offset=0.1, font=2)
          lines(c(max(MAT.change), max(MAT.change)+position[i]*0.6), rep(line[length(line)],2), lty=2, lwd=0.6)
        } else {
          text(seq(0,max(MAT.change), 0.01)[which(line==max(line))]+position[i]*0.2, max(line), spp, pos=3, cex=0.8, offset=0.1, font=2)
        }
      }
    }
    rect(0,0,max(MAT.change),ylim[2]*1.1, col=NA, border=T)
    box()
    
    # boxplot for focal period
    par(xpd=T)
    rcp.focal <- "rcp45"
    for(proj.year.focal in proj.years){
      x <- c(0,MAT.change)
      x.focal <- MAT.change[which(seq.rcp==rcp.focal & seq.proj.year==proj.year.focal)]
      position <- ylim[1] - diff(ylim)/40 - diff(ylim)/37.5*which(proj.years==proj.year.focal)
      boxplot(x.focal, add=T, horizontal=TRUE, axes=FALSE, range=0, at=position, boxwex = diff(ylim)/30)
      text(max(x.focal), position, paste(rcp.name[which(rcps==rcp.focal)], ", ", proj.year.name[which(proj.years==proj.year.focal)], sep=""), pos=4, cex=0.8)
    }
    par(xpd=F)
    
    mtext(if(transform==T) "(A)" else "(B)", side=3, line=-1.5, adj=0.25, cex=1, font=2)
    mtext("Native", side=3, line=-1.5, adj=0.05, cex=1, font=2)
    mtext("Exotic", side=3, line=-1.5, adj=0.95, cex=1, font=2)
    
  }
  dev.off() 
  print(edatope)
}




############################
## Three panel plot of species trends relative to MAT change, by edatope
############################
native <- spps.lookup$Native[match(spps,spps.lookup$TreeCode)]
spps.native <- c("Ra", "Pl", "Pj", "Fd", "Cw", "Ba", "Sx", "Bl", "Bg", "Yc", "Pa", "Hm", "Lw", "La", "Lt", "Hw", "Py", "Dr", "Ep", "At", "Acb", "Pw", "Ss", "Sb", "Qg", "Act", "Mb")

  
  png(filename=paste("Results\\ManuscriptFigures\\Manu_Spaghetti\\CCISS_manu_SppSpaghetti.png",sep="."), type="cairo", units="in", width=6.5, height=5, pointsize=10, res=400)
  mat <- matrix(c(1,2,3,4, 6, 5,5,5),2, byrow=T)   #define the plotting order
  layout(mat, widths=c(0.25,1,1,1), heights=c(1, 0.05))   #set up the multipanel plot
  
  par(mar=c(0,0,0,0))
  plot(1, type="n", axes=F, xlab="", ylab="")  
  text(0.75,1,"Suitable area (sq. km)", srt=90, cex=1.2)
  
  for(edatope in edatopes){
    for(spp in spps){
    x <- c(0,MAT.change)
    y <- c(get(paste("SuitCells.ref", spp, edatope, sep=".")),get(paste("SuitCells", spp, edatope, sep=".")))*4 #times 4km^2 because they are 2km grid cells. 
    l <- loess(y[order(x)]~x[order(x)])
    assign(paste("line",spp, sep="."), predict(l, seq(0,max(MAT.change), 0.01)))
  }
  
      transform=T
  # for(transform in c(T, F)){
    
    par(mar=c(3.1,0,0,0.2), mgp=c(4, 0.2, 0))
    ylim=if(transform==T) c(2.4,5.9) else c(0,47000)
    plot(0, xlim=c(-2.5,9.5), ylim=ylim, yaxs="i", xaxs="i", col="white", xaxt="n", yaxt="n", 
         xlab="", 
         ylab="")
    # if(edatope==edatopes[2]){
    #   par(xpd=T)
    #   title(xlab=list(bquote(BC~mean~temperature~change~relative~to~"1961-90"~"("*degree*C*")"), cex=1.2))
    #   par(xpd=F)
    # }
    axis(1, at=0:7, labels = 0:7, tck=0)
    if(edatope==edatopes[1]){
    par(mgp=c(1, 0.2, 0))
    par(xpd=T)
    axis(2, lty=0, at=if(transform==T) log10(c(1000, 5000, 10000, 50000, 100000, 500000)) else seq(0,50000,10000), labels = if(transform==T) format(c(1000, 5000, 10000, 50000, 100000, 500000), scientific = FALSE, big.mark=",") else format(seq(0,50000,10000), scientific = FALSE, big.mark=","), tck=0, las=2)
    par(xpd=F)
    }
    # rect(-9,0,0, 60000, col="lightgray", border=F)
    # rect(max(MAT.change),0,9, ylim[2]*1.1, col="lightgray", border=F)
    
    suit.exotic.final <- vector()
    for(spp in spps[-which(spps%in%spps.native)]){
      line <- get(paste("line",spp, sep="."))
      if(transform==T) line[line<1] <- 1
      if(transform==T) line <- log10(line)
      suit.exotic.final[which(spps[-which(spps%in%spps.native)]==spp)] <- line[length(line)]
    }
    
    suit.native.initial <- vector()
    for(spp in spps[which(spps%in%spps.native)]){
      suit.native.initial[which(spps[which(spps%in%spps.native)]==spp)] <- get(paste("SuitCells.ref", spp, edatope, sep="."))
    }
    
    spplist <- spps[which(spps%in%spps.native)][order(suit.native.initial)]
    colors = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)][-1]
    colors = colors[-grep("yellow", colors)]
    set.seed(5)
    ColScheme <- c(brewer.pal(n=12, "Paired")[-11],sample(colors,length(spps)-11))
    for(spp in spplist){
      i <- which(spplist==spp)
      line <- get(paste("line",spp, sep="."))
      if(transform==T) line[line<1] <- 1
      if(transform==T) line <- log10(line)
      if(line[1]> if(transform==T) ylim[1] else 100){
        lines(seq(0,max(MAT.change), 0.01), line, col=ColScheme[i], lwd=2)
        position <- rep(0:3, times=100)
        text(0-position[i]*0.55, line[1], spp, pos=2, col=ColScheme[i], font=2, cex=0.9, offset=0.1)
        lines(c(0-position[i]*0.55,0), rep(line[1],2), col=ColScheme[i], lty=2)
      }
    }
    
    spplist <- spps[-which(spps%in%spps.native)][order(suit.exotic.final)]
    for(spp in spplist){
      i <- which(spplist==spp)
      line <- get(paste("line",spp, sep="."))
      if(transform==T) line[line<1] <- 1
      if(transform==T) line <- log10(line)
      if(max(line)> if(transform==T) ylim[1] else 100){
        lines(seq(0,max(MAT.change), 0.01), line)
        position <- rep(0:3, times=100)
        if(which.max(line)>(length(line)-100)){
          text(max(MAT.change)+position[i]*0.55, line[length(line)], spp, pos=4, cex=0.9, offset=0.1, font=2)
          lines(c(max(MAT.change), max(MAT.change)+position[i]*0.55), rep(line[length(line)],2), lty=2, lwd=0.6)
        } else {
          text(seq(0,max(MAT.change), 0.01)[which(line==max(line))]+position[i]*0.2, max(line), spp, pos=3, cex=0.9, offset=0.1, font=2)
        }
      }
    }
    rect(0,0,max(MAT.change),ylim[2]*1.1, col=NA, border=T)
    box()
    
    
    # boxplot for focal period
    if(edatope==edatopes[2]){
    par(xpd=T)
    rcp.focal <- "rcp45"
    for(proj.year.focal in proj.years){
      x <- c(0,MAT.change)
      x.focal <- MAT.change[which(seq.rcp==rcp.focal & seq.proj.year==proj.year.focal)]
      position <- ylim[1] - diff(ylim)/50 - diff(ylim)/60*which(proj.years==proj.year.focal)
      boxplot(x.focal, add=T, horizontal=TRUE, axes=FALSE, range=0, at=position, boxwex = diff(ylim)/50)
      text(max(x.focal), position, paste(rcp.name[which(rcps==rcp.focal)], ", ", proj.year.name[which(proj.years==proj.year.focal)], sep=""), pos=4, cex=0.9)
    }
    par(xpd=F)
  }
    
    mtext(paste("(", LETTERS[which(edatopes==edatope)],") ", edatope, " sites", sep=""), side=3, line=-1.5, adj=0.325, cex=0.8, font=2)
    mtext("Native", side=3, line=-1.5, adj=0.025, cex=0.8, font=2)
    mtext("Exotic", side=3, line=-1.5, adj=0.975, cex=0.8, font=2)
    
  # }
  print(edatope)
  
  
  }
  par(mar=c(0,0,0,0))
  plot(1, type="n", axes=F, xlab="", ylab="")  
  text(1,1,bquote(BC~mean~temperature~change~relative~to~"1961-90"~"("*degree*C*")"), srt=0, cex=1.2)

    dev.off() 
  

