
##======================================================================================
## CCISS Publication Scripts
## Step 4c - Figure - Spaghetti plots of spp suitable area relative to MAT change
##======================================================================================

# Colin Mahony
# c_mahony@alumni.ubc.ca
# 778-288-4008
# July 21, 2019



#===============================================================================
# Set analysis Parameters
#===============================================================================

source("./_CCISS_Packages.R") ## packages required
source("./_CCISS_Functions.R") ## common functions
source("./_CCISS_Parameters.R") ## settings used through all scripts

## non-THLB BGCs for exclusion from results
points <- read.csv(paste("inputs\\",grid,".csv", sep=""))
BGC <- points$ID2
BGC <- gsub(" ","",BGC)
BGCs_notin_THLB <- BGCs_notin_THLB$BGC[which(BGCs_notin_THLB$Exlude=="x")]
exclude <- which(BGC%in%BGCs_notin_THLB)

#===============================================================================
# Import suitability tables
#===============================================================================
S1 <- treesuit
S1 <- unique(S1)[,-5]
dim(S1)
S1 <- unique(S1)
dim(S1)

# select a subset of species to run the analysis on
spps <- unique(S1$Spp)
spps.candidate <- spps.lookup$TreeCode[-which(spps.lookup$Exclude=="x")]
spps <- spps[which(spps%in%spps.candidate)]

native <- spps.lookup$Native[match(spps,spps.lookup$TreeCode)]
spps.native <- c("Ra", "Pl", "Pj", "Fd", "Cw", "Ba", "Sx", "Bl", "Bg", "Yc", "Pa", "Hm", "Lw", "La", "Lt", "Hw", "Py", "Dr", "Ep", "At", "Acb", "Pw", "Ss", "Sb", "Qg", "Act", "Mb")




#===============================================================================
# calculate mean MAT change for each model prediction
#===============================================================================

fplot=paste("inputs\\", grid, "_Normal_1961_1990MSY.csv", sep="")
Y0 <- fread(fplot, select = "MAT", stringsAsFactors = FALSE, data.table = FALSE) #fread is faster than read.csv
MAT.ref <- Y0$MAT
MAT.mean.ref <- mean(MAT.ref, na.rm=T)

for(hist.year in hist.years){
  Y0 <- fread(paste("inputs\\", grid, "_", hist.year, "_BioVars.csv", sep=""), select = "MAT", stringsAsFactors = FALSE, data.table = FALSE) #fread is faster than read.csv
  assign(paste("MAT", hist.year, sep="."), Y0$MAT)
  assign(paste("MAT.change", hist.year, sep="."), mean(Y0$MAT, na.rm=T)-MAT.mean.ref)
  print(hist.year)
}

for(GCM in GCMs){
  Y1 <- fread(paste("inputs\\", grid, "_", GCM, ".csv", sep=""), select = c("Year", "MAT"), stringsAsFactors = FALSE, data.table = FALSE)
  ## assign single vectors to RCPs and proj.years
  Ystr <- strsplit(Y1[,1], "_")
  Y4 <- matrix(unlist(Ystr), ncol=3, byrow=TRUE)
  Y4[,3] <- gsub(".gcm","",Y4[,3])
  for(rcp in rcps){
    for(proj.year in proj.years){
      assign(paste("MAT", GCM, rcp, proj.year, sep="."), Y1$MAT[which(Y4[,2]==rcp & Y4[,3]==proj.year)])
    }
    # print(rcp)
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
    Suit <- read.csv(paste("outputs\\Suit.ref", grid, spp, edatope, "csv", sep="."))[-exclude,1]
    Suit[is.na(Suit)] <- 5
    Suit <- 1-(Suit-1)/4
    assign(paste("SuitCells.ref", spp, edatope, sep="."), sum(Suit))
    assign(paste("SppCells.ref", spp, edatope, sep="."), sum(Suit>0))
    # print(edatope)
  }
  print(spp)
}

n.cells <- length(Suit)

for(spp in spps){
  for(edatope in edatopes){
    for(rcp in rcps){
      for(proj.year in proj.years){
        
        SuitCells <- rep(NA, length(GCMs))
        SppCells <- rep(NA, length(GCMs))
        for(GCM in GCMs){
          
          Suit.proj <- read.csv(paste("outputs\\Suit", grid, GCM, rcp, proj.year, spp, edatope, "csv", sep="."))[-exclude,1]
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
rcp.focal <- "rcp45"
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
## Three panel plot of species trends relative to MAT change, by edatope
############################


png(filename=paste("results\\Manu_Spaghetti\\CCISS_manu_SppSpaghetti.png",sep="."), type="cairo", units="in", width=6.5, height=5, pointsize=10, res=400)
# pdf(file=paste("results\\CCISS.Fig5.SppSpaghetti","pdf",sep="."), width=6.5, height=5, pointsize=10)
mat <- matrix(c(1,2,3,4, 6, 5,5,5),2, byrow=T)   #define the plotting order
layout(mat, widths=c(0.225,1,1,1), heights=c(1, 0.04))   #set up the multipanel plot

par(mar=c(0,0,0,0))
plot(1, type="n", axes=F, xlab="", ylab="")  
text(0.75,1,"Feasible range (% of edatope area in BC)", srt=90, cex=1.2)

for(edatope in edatopes){
  for(spp in spps){
    x <- c(0,MAT.change)
    # y <- c(get(paste("SuitCells.ref", spp, edatope, sep=".")),get(paste("SuitCells", spp, edatope, sep=".")))*4 #times 4km^2 because they are 2km grid cells. 
    y <- c(get(paste("SuitCells.ref", spp, edatope, sep=".")),get(paste("SuitCells", spp, edatope, sep=".")))/n.cells #divided by number of cells in BC to get proportion of edatope area
    l <- loess(y[order(x)]~x[order(x)])
    assign(paste("line",spp, sep="."), predict(l, seq(0,max(MAT.change), 0.01)))
  }
  
  transform=T
  # for(transform in c(T, F)){
  
  par(mar=c(3.85,0,0,0.2), mgp=c(4, 0.2, 0))
  ylim=if(transform==T) c(-2.65,0.01) else c(0,1)
  plot(0, xlim=c(-3.5,8.25), ylim=ylim, yaxs="i", xaxs="i", col="white", xaxt="n", yaxt="n", 
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
    y.labels <- c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5)
    axis(2, lty=0, at=if(transform==T) log10(y.labels) else seq(0,1,0.1), labels = if(transform==T) paste(y.labels*100, "%", sep="") else format(seq(0,1,0.1), scientific = FALSE, big.mark=","), tck=0, las=2)
    par(xpd=F)
  }
  # rect(-9,0,0, 60000, col="lightgray", border=F)
  # rect(max(MAT.change),0,9, ylim[2]*1.1, col="lightgray", border=F)
  
  suit.exotic.final <- vector()
  for(spp in spps[-which(spps%in%spps.native)]){
    line <- get(paste("line",spp, sep="."))
    if(transform==T) line[line<10^(-5)] <- 10^(-5)
    if(transform==T) line <- log10(line)
    suit.exotic.final[which(spps[-which(spps%in%spps.native)]==spp)] <- line[length(line)]
  }
  
  suit.native.initial <- vector()
  for(spp in spps[which(spps%in%spps.native)]){
    suit.native.initial[which(spps[which(spps%in%spps.native)]==spp)] <- get(paste("SuitCells.ref", spp, edatope, sep="."))
  }
  
  spplist <- spps[which(spps%in%spps.native)][order(suit.native.initial)]
  
  # #Color scheme for individual species
  # colors = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)][-1]
  # colors = colors[-grep("yellow", colors)]
  # set.seed(5)
  # ColScheme <- c(brewer.pal(n=12, "Paired")[-11],sample(colors,length(spps)-11))
  
  # #Color scheme for species groups
  boreal <- c("Pl", "Sx", "Sb", "At", "Ep", "Pj", "Acb")
  temperate <- c("Fd", "Lw", "Pw", "Py", "Bg", "Act")
  mesothermal <- c("Hw", "Cw", "Ba", "Ss", "Dr", "Mb")
  subalpine <- c("Hm", "Yc", "Bl", "Ba")
  ColScheme <- rep(NA, length(spplist))
  ColScheme[which(spplist%in%boreal)] <- as.character(BGCcolors.BC$HEX[which(BGCcolors.BC$zone=="SBS")])
  ColScheme[which(spplist%in%temperate)] <- as.character(BGCcolors.BC$HEX[which(BGCcolors.BC$zone=="IDF")])
  ColScheme[which(spplist%in%mesothermal)] <- as.character(BGCcolors.BC$HEX[which(BGCcolors.BC$zone=="CWH")])
  ColScheme[which(spplist%in%subalpine)] <- as.character(BGCcolors.BC$HEX[which(BGCcolors.BC$zone=="MS")])
  
  if(edatope==edatopes[2]){
    text(-3.2, ylim[1]+0.02, "Boreal species", cex=1.1, srt=90, font=2, pos=4, col=unique(ColScheme[which(spplist%in%boreal)]))
    text(-2.4, ylim[1]+0.02, "Temperate species", cex=1.1, srt=90, font=2, pos=4, col=unique(ColScheme[which(spplist%in%temperate)]))
    text(-1.6, ylim[1]+0.02, "Mesothermal species", cex=1.1, srt=90, font=2, pos=4, col=unique(ColScheme[which(spplist%in%mesothermal)]))
    text(-0.8, ylim[1]+0.02, "Subalpine species", cex=1.1, srt=90, font=2, pos=4, col=unique(ColScheme[which(spplist%in%subalpine)]))
  }
  
  for(spp in spplist){
    i <- which(spplist==spp)
    line <- get(paste("line",spp, sep="."))
    if(transform==T) line[line<10^(-5)] <- 10^(-5)
    if(transform==T) line <- log10(line)
    if(line[1]> if(transform==T) ylim[1] else 100){
      lines(seq(0,max(MAT.change), 0.01), line, col=ColScheme[i], lwd=2)
      position <- rep(0:3, times=100)
      text(0-position[i]*0.8, line[1], spp, pos=2, col=ColScheme[i], font=2, cex=0.9, offset=0.1)
      lines(c(0-position[i]*0.8,0), rep(line[1],2), col=ColScheme[i], lty=2)
    }
  }
  
  spplist <- spps[-which(spps%in%spps.native)][order(suit.exotic.final)]
  for(spp in spplist){
    i <- which(spplist==spp)
    line <- get(paste("line",spp, sep="."))
    if(transform==T) line[line<10^(-5)] <- 10^(-5)
    if(transform==T) line <- log10(line)
    if(max(line)> if(transform==T) ylim[1] else 100){
      lines(seq(0,max(MAT.change), 0.01), line)
      position <- rep(0, times=100)
      if(which.max(line)>(length(line)-100)){
        text(max(MAT.change)+position[i]*0.55, line[length(line)], spp, pos=4, cex=0.9, offset=0.1, font=2)
        lines(c(max(MAT.change), max(MAT.change)+position[i]*0.55), rep(line[length(line)],2), lty=2, lwd=0.6)
      } else {
        text(seq(0,max(MAT.change), 0.01)[which(line==max(line))]+position[i]*0.2, max(line), spp, pos=3, cex=0.9, offset=0.1, font=2)
      }
    }
  }
  rect(0,-10,max(MAT.change),ylim[2]*1.1, col=NA, border=T)
  box()
  
  
  # boxplots
  par(xpd=T)
  if(edatope==edatopes[2]){
    for(i in 1:4){
      rcp.focal <- rcps[c(1,1,1,2)][i]
      proj.year.focal <- proj.years[c(1,2,3,3)][i]
      x <- MAT.change
      x.focal <- MAT.change[which(seq.rcp[-1]==rcp.focal & seq.proj.year[-1]==proj.year.focal)]
      position <- ylim[1] - diff(ylim)/40 - diff(ylim)/60*i
      boxplot(x.focal, add=T, col=c("dodgerblue", "red")[which(rcps==rcp.focal)], horizontal=TRUE, axes=FALSE, range=0, at=position, boxwex = diff(ylim)/50)
      text(if(rcp.focal=="rcp85" & proj.year.focal==2085) min(x.focal) else max(x.focal), position, paste(rcp.name[which(rcps==rcp.focal)], ", ", proj.year.name[which(proj.years==proj.year.focal)], sep=""), pos=if(rcp.focal=="rcp85" & proj.year.focal==2085) 2 else 4, cex=0.9)
      
    }
  }
  
  mtext(paste("(", letters[which(edatopes==edatope)],") ", edatope, " edatope", sep=""), side=3, line=-1.5, adj=0.5, cex=0.8, font=2)
  mtext("Native", side=3, line=-1.5, adj=0.025, cex=0.8, font=2)
  mtext("Non-native", side=4, line=-1.5, adj=0.975, cex=0.8, font=2)
  
  par(xpd=F)
  
  # }
  print(edatope)
  
  
}
par(mar=c(0,0,0,0))
plot(1, type="n", axes=F, xlab="", ylab="")  
text(1,1,bquote(BC~mean~temperature~change~relative~to~"1961-90"~"("*degree*C*")"), srt=0, cex=1.2)

dev.off()



# ############################
# ## two-panel plot of species trends relative to MAT change, log and raw scaled
# ############################
# 
# edatope="C4"
# for(edatope in edatopes){
# 
#   png(filename=paste("results\\CCISS_manu_SppSpaghetti", edatope, "png",sep="."), type="cairo", units="in", width=6.5, height=5, pointsize=8, res=400)
# 
#   for(spp in spps){
#     x <- c(0,MAT.change)
#     y <- c(get(paste("SuitCells.ref", spp, edatope, sep=".")),get(paste("SuitCells", spp, edatope, sep=".")))*4 #times 4km^2 because they are 2km grid cells.
#     l <- loess(y[order(x)]~x[order(x)])
#     assign(paste("line",spp, sep="."), predict(l, seq(0,max(MAT.change), 0.01)))
#   }
# 
#   par(mfrow=c(1,2))
#   for(transform in c(T, F)){
# 
#     par(mar=c(5,4,0,0), mgp=c(4, 0.2, 0))
#     ylim=if(transform==T) c(2,6.1) else c(0,47000)
#     plot(0, xlim=c(-2.5,9.5), ylim=ylim, yaxs="i", xaxs="i", col="white", xaxt="n", yaxt="n",
#          xlab=bquote(BC~mean~temperature~change~relative~to~"1961-90"~"("*degree*C*")"),
#          ylab="")
#     axis(1, at=0:7, labels = 0:7, tck=0)
#     par(mgp=c(3, 0.2, 0))
#     if(transform==T) title(ylab="Feasible area (sq. km)")
#     axis(2, at=if(transform==T) 0:7 else seq(0,50000,10000), labels = if(transform==T) format(10^(0:7), scientific = FALSE, big.mark=",") else format(seq(0,50000,10000), scientific = FALSE, big.mark=","), tck=0, las=2)
#     # rect(-9,0,0, 60000, col="lightgray", border=F)
#     # rect(max(MAT.change),0,9, ylim[2]*1.1, col="lightgray", border=F)
# 
#     suit.exotic.final <- vector()
#     for(spp in spps[-which(spps%in%spps.native)]){
#       line <- get(paste("line",spp, sep="."))
#       if(transform==T) line[line<1] <- 1
#       if(transform==T) line <- log10(line)
#       suit.exotic.final[which(spps[-which(spps%in%spps.native)]==spp)] <- line[length(line)]
#     }
# 
#     suit.native.initial <- vector()
#     for(spp in spps[which(spps%in%spps.native)]){
#       suit.native.initial[which(spps[which(spps%in%spps.native)]==spp)] <- get(paste("SuitCells.ref", spp, edatope, sep="."))
#     }
# 
#     spplist <- spps[which(spps%in%spps.native)][order(suit.native.initial)]
#     colors = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)][-1]
#     set.seed(2)
#     ColScheme <- c(brewer.pal(n=12, "Paired"),sample(colors,length(spps)-12))
#     for(spp in spplist){
#       i <- which(spplist==spp)
#       line <- get(paste("line",spp, sep="."))
#       if(transform==T) line[line<1] <- 1
#       if(transform==T) line <- log10(line)
#       if(line[1]> if(transform==T) ylim[1] else 100){
#         lines(seq(0,max(MAT.change), 0.01), line, col=ColScheme[i], lwd=2)
#         position <- rep(0:3, times=100)
#         text(0-position[i]*0.6, line[1], spp, pos=2, col=ColScheme[i], font=2, cex=0.75, offset=0.1)
#         lines(c(0-position[i]*0.6,0), rep(line[1],2), col=ColScheme[i], lty=2)
#       }
#     }
# 
#     spplist <- spps[-which(spps%in%spps.native)][order(suit.exotic.final)]
#     for(spp in spplist){
#       i <- which(spplist==spp)
#       line <- get(paste("line",spp, sep="."))
#       if(transform==T) line[line<1] <- 1
#       if(transform==T) line <- log10(line)
#       if(max(line)> if(transform==T) ylim[1] else 100){
#         lines(seq(0,max(MAT.change), 0.01), line)
#         position <- rep(0:3, times=100)
#         if(which.max(line)>(length(line)-50)){
#           text(max(MAT.change)+position[i]*0.6, line[length(line)], spp, pos=4, cex=0.75, offset=0.1, font=2)
#           lines(c(max(MAT.change), max(MAT.change)+position[i]*0.6), rep(line[length(line)],2), lty=2, lwd=0.6)
#         } else {
#           text(seq(0,max(MAT.change), 0.01)[which(line==max(line))]+position[i]*0.2, max(line), spp, pos=3, cex=0.8, offset=0.1, font=2)
#         }
#       }
#     }
#     rect(0,0,max(MAT.change),ylim[2]*1.1, col=NA, border=T)
#     box()
# 
#     # boxplot for focal period
#     par(xpd=T)
#     rcp.focal <- "rcp45"
#     for(proj.year.focal in proj.years){
#       x <- c(0,MAT.change)
#       x.focal <- MAT.change[which(seq.rcp==rcp.focal & seq.proj.year==proj.year.focal)]
#       position <- ylim[1] - diff(ylim)/40 - diff(ylim)/37.5*which(proj.years==proj.year.focal)
#       boxplot(x.focal, add=T, horizontal=TRUE, axes=FALSE, range=0, at=position, boxwex = diff(ylim)/30)
#       text(max(x.focal), position, paste(rcp.name[which(rcps==rcp.focal)], ", ", proj.year.name[which(proj.years==proj.year.focal)], sep=""), pos=4, cex=0.8)
#     }
#     par(xpd=F)
# 
#     mtext(if(transform==T) "(a)" else "(b)", side=3, line=-1.5, adj=0.25, cex=1, font=2)
#     mtext("Native", side=3, line=-1.5, adj=0.05, cex=1, font=2)
#     mtext("Exotic", side=3, line=-1.5, adj=0.95, cex=1, font=2)
# 
#   }
#   dev.off()
#   print(edatope)
# }

############################
## separate out species groups
############################


edatope="C4"
for(edatope in edatopes){
  for(spp in spps){
    x <- c(0,MAT.change)
    # y <- c(get(paste("SuitCells.ref", spp, edatope, sep=".")),get(paste("SuitCells", spp, edatope, sep=".")))*4 #times 4km^2 because they are 2km grid cells. 
    y <- c(get(paste("SuitCells.ref", spp, edatope, sep=".")),get(paste("SuitCells", spp, edatope, sep=".")))/n.cells #divided by number of cells in BC to get proportion of edatope area
    l <- loess(y[order(x)]~x[order(x)])
    assign(paste("line",spp, sep="."), predict(l, seq(0,max(MAT.change), 0.01)))
  }
  
  transform=T
  # for(transform in c(T, F)){
  
  suit.exotic.final <- vector()
  for(spp in spps[-which(spps%in%spps.native)]){
    line <- get(paste("line",spp, sep="."))
    if(transform==T) line[line<10^(-5)] <- 10^(-5)
    if(transform==T) line <- log10(line)
    suit.exotic.final[which(spps[-which(spps%in%spps.native)]==spp)] <- line[length(line)]
  }
  
  suit.native.initial <- vector()
  for(spp in spps[which(spps%in%spps.native)]){
    suit.native.initial[which(spps[which(spps%in%spps.native)]==spp)] <- get(paste("SuitCells.ref", spp, edatope, sep="."))
  }
  
  spplist <- spps[which(spps%in%spps.native)][order(suit.native.initial)]
  
  # #Color scheme for individual species
  # colors = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)][-1]
  # colors = colors[-grep("yellow", colors)]
  # set.seed(5)
  # ColScheme <- c(brewer.pal(n=12, "Paired")[-11],sample(colors,length(spps)-11))
  
  # #Color scheme for species groups
  spp.types <- c("boreal", "temperate", "mesothermal", "subalpine", "exotic")
  boreal <- c("Pl", "Sx", "Sb", "At", "Ep", "Pj", "Acb")
  temperate <- c("Fd", "Lw", "Pw", "Py", "Bg", "Act")
  mesothermal <- c("Hw", "Cw", "Ss", "Dr", "Mb")
  subalpine <- c("Hm", "Yc", "Bl", "Ba")
  ColScheme <- rep(NA, length(spplist))
  ColScheme[which(spplist%in%boreal)] <- as.character(BGCcolors.BC$HEX[which(BGCcolors.BC$zone=="SBS")])
  ColScheme[which(spplist%in%temperate)] <- as.character(BGCcolors.BC$HEX[which(BGCcolors.BC$zone=="ICH")])
  ColScheme[which(spplist%in%mesothermal)] <- as.character(BGCcolors.BC$HEX[which(BGCcolors.BC$zone=="CWH")])
  ColScheme[which(spplist%in%subalpine)] <- as.character(BGCcolors.BC$HEX[which(BGCcolors.BC$zone=="MS")])
  
  # if(edatope==edatopes[2]){
  #   text(-3.2, ylim[1]+0.02, "Boreal species", cex=1.1, srt=90, font=2, pos=4, col=unique(ColScheme[which(spplist%in%boreal)]))
  #   text(-2.4, ylim[1]+0.02, "Temperate species", cex=1.1, srt=90, font=2, pos=4, col=unique(ColScheme[which(spplist%in%temperate)]))
  #   text(-1.6, ylim[1]+0.02, "Mesothermal species", cex=1.1, srt=90, font=2, pos=4, col=unique(ColScheme[which(spplist%in%mesothermal)]))
  #   text(-0.8, ylim[1]+0.02, "Subalpine species", cex=1.1, srt=90, font=2, pos=4, col=unique(ColScheme[which(spplist%in%subalpine)]))
  # }
  
  spp.type.names <- c("Boreal", "Temperate", "Mesothermal", "Subalpine", "Exotic")
  
  for(k in 1:length(spp.types)){
    spp.type <- spp.types[k]
    
    png(filename=paste("results\\Manu_Spaghetti\\CCISS_SppSpaghetti.spptypes", spp.type, edatope, "png",sep="."), type="cairo", units="in", width=3.5, height=5, pointsize=10, res=400)
    # pdf(file=paste("results\\CCISS.Fig5.SppSpaghetti","pdf",sep="."), width=6.5, height=5, pointsize=10)
    mat <- matrix(c(1,2,4,3),2, byrow=T)   #define the plotting order
    layout(mat, widths=c(0.225,1), heights=c(1, 0.04))   #set up the multipanel plot
    
    par(mar=c(0,0,0,0))
    plot(1, type="n", axes=F, xlab="", ylab="")  
    text(0.75,1,paste("Feasible range (% of", if(edatope=="C4") "zonal" else "edatope","site area in BC)"), srt=90, cex=1.2)
    par(mar=c(3.85,0,2,0.2), mgp=c(4, 0.2, 0))
  ylim=if(transform==T) c(-2.15,0.01) else c(0,1)
  plot(0, xlim=c(-3.5,8.25), ylim=ylim, yaxs="i", xaxs="i", col="white", xaxt="n", yaxt="n", 
       xlab="", 
       ylab="")
  # if(edatope==edatopes[2]){
  #   par(xpd=T)
  #   title(xlab=list(bquote(BC~mean~temperature~change~relative~to~"1961-90"~"("*degree*C*")"), cex=1.2))
  #   par(xpd=F)
  # }
  axis(1, at=0:7, labels = 0:7, tck=0)
  par(mgp=c(1, 0.2, 0))
  par(xpd=T)
  y.labels <- c(0.01, 0.025, 0.05, 0.1, 0.25, 0.5)
  axis(2, lty=0, at=if(transform==T) log10(y.labels) else seq(0,1,0.1), labels = if(transform==T) paste(y.labels*100, "%", sep="") else format(seq(0,1,0.1), scientific = FALSE, big.mark=","), tck=0, las=2)
  par(xpd=F)
  # rect(-9,0,0, 60000, col="lightgray", border=F)
  # rect(max(MAT.change),0,9, ylim[2]*1.1, col="lightgray", border=F)
  
  
  if(k<5){
    for(spp in get(spp.types[k])){
    i <- which(spplist==spp)
    line <- get(paste("line",spp, sep="."))
    if(transform==T) line[line<10^(-5)] <- 10^(-5)
    if(transform==T) line <- log10(line)
    if(line[1]> if(transform==T) ylim[1] else 100){
      lines(seq(0,4, 0.01), line[which(seq(0,max(MAT.change), 0.01)<=4)], col=ColScheme[i], lwd=2)
      lines(seq(4,max(MAT.change), 0.01), line[which(seq(0,max(MAT.change), 0.01)>=4)], col=ColScheme[i], lwd=2, lty=2)
      position <- rep(0:3, times=100)
      text(0-position[i]*0.8, line[1], spp, pos=2, col=ColScheme[i], font=2, cex=0.9, offset=0.1)
      lines(c(0-position[i]*0.8,0), rep(line[1],2), col=ColScheme[i], lty=2)
      text(max(MAT.change), line[length(line)], spp, pos=4, col=ColScheme[i], font=2, cex=0.9, offset=0.1)
    }
  }
  }
  
  if(k==5){
  spplist <- spps[-which(spps%in%spps.native)][order(suit.exotic.final)]
  for(spp in spplist){
    i <- which(spplist==spp)
    line <- get(paste("line",spp, sep="."))
    if(transform==T) line[line<10^(-5)] <- 10^(-5)
    if(transform==T) line <- log10(line)
    if(max(line)> if(transform==T) ylim[1] else 100){
      lines(seq(0,4, 0.01), line[which(seq(0,max(MAT.change), 0.01)<=4)], lwd=2)
      lines(seq(4,max(MAT.change), 0.01), line[which(seq(0,max(MAT.change), 0.01)>=4)], lwd=2, lty=2)
      position <- rep(0, times=100)
      if(which.max(line)>(length(line)-100)){
        text(max(MAT.change)+position[i]*0.55, line[length(line)], spp, pos=4, cex=0.9, offset=0.1, font=2)
        lines(c(max(MAT.change), max(MAT.change)+position[i]*0.55), rep(line[length(line)],2), lty=2, lwd=0.6)
      } else {
        text(seq(0,max(MAT.change), 0.01)[which(line==max(line))]+position[i]*0.2, max(line), spp, pos=3, cex=0.9, offset=0.1, font=2)
      }
    }
  }
  }
  
  rect(0,-10,max(MAT.change),ylim[2]*1.1, col=NA, border=T)
  box()
  
  
  # boxplots
  par(xpd=T)
    for(i in 1:4){
      rcp.focal <- rcps[c(1,1,1,2)][i]
      proj.year.focal <- proj.years[c(1,2,3,3)][i]
      x <- MAT.change
      x.focal <- MAT.change[which(seq.rcp[-1]==rcp.focal & seq.proj.year[-1]==proj.year.focal)]
      position <- ylim[1] - diff(ylim)/35 - diff(ylim)/45*i
      boxplot(x.focal, add=T, col=c("dodgerblue", "red")[which(rcps==rcp.focal)], horizontal=TRUE, axes=FALSE, range=0, at=position, boxwex = diff(ylim)/50)
      text(if(rcp.focal=="rcp85" & proj.year.focal==2085) min(x.focal) else max(x.focal), position, paste(rcp.name[which(rcps==rcp.focal)], ", ", proj.year.name[which(proj.years==proj.year.focal)], sep=""), pos=if(rcp.focal=="rcp85" & proj.year.focal==2085) 2 else 4, cex=0.9)
      
    }

  # mtext(paste("(", letters[which(edatopes==edatope)],") ", edatope, " edatope", sep=""), side=3, line=-1.5, adj=0.5, cex=0.8, font=2)
  mtext(paste(spp.type.names[k], "species"), side=3, line=0.25, adj=0.025, cex=1.1, font=2)
  # mtext("Non-native", side=4, line=-1.5, adj=0.975, cex=0.8, font=2)
  
  par(xpd=F)
  
  # }
  par(mar=c(0,0,0,0))
  plot(1, type="n", axes=F, xlab="", ylab="")  
  text(1,1,bquote(BC~temperature~change~relative~to~"1961-90"~"("*degree*C*")"), srt=0, cex=1)
  
  
  dev.off()
  print(spp.types[k])
}
print(edatope)
}
