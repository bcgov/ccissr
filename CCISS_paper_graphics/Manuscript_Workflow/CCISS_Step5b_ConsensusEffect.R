
##======================================================================================
## CCISS Publication Scripts
## Step 5b - analysis of whether there is a bias in the "consensus" projection [there isn't]
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

rcp.focal="rcp45"
proj.year.focal=2025

#BGC zone color scheme
BGCcolors$colour <- as.character(BGCcolors$colour)
BGCcolors$colour[match(BGCcolors.BC$zone, BGCcolors$classification)] <- as.character(BGCcolors.BC$HEX)
ColScheme <- factor(BGCcolors$colour, levels=BGCcolors$colour)
zones <- factor(BGCcolors$classification, levels=BGCcolors$classification)

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

#===============================================================================
# BGC change for each model prediction
#===============================================================================

## mapped BGC
points <- read.csv(paste("inputs\\",grid,".csv", sep=""))
BGC <- points$ID2
BGC <- gsub(" ","",BGC)  
zone <- rep(NA, length(BGC))
for(i in BGCcolors.BC$zone){ zone[grep(i,BGC)] <- i }

## reference period BGC
BGC.pred.ref <- as.character(read.csv(paste("outputs\\BGC.pred", grid, "ref", model, "csv", sep="."), header=F)[,1])
zone.pred.ref <- rep(NA, length(BGC))
for(i in BGCcolors.BC$zone){ zone.pred.ref[grep(i,BGC.pred.ref)] <- i }

# Historical BGC
for(hist.year in hist.years){
  BGC.pred <- as.character(read.csv(paste("outputs\\BGC.pred", grid,hist.year, model,"csv", sep="."), header=F)[,1])
  assign(paste("BGC.pred", hist.year, sep="."), BGC.pred) #bgc projection
  
  #did the BGC change?
  assign(paste("BGC.change", hist.year, sep="."), BGC.pred!=BGC) #did the BGC unit change (mapped baseline)
  assign(paste("BGC.change.pred", hist.year, sep="."), BGC.pred!=BGC.pred.ref) #did the BGC unit change (predicted baseline)
  
  #did the zone change? 
  zone.pred <- rep(NA, length(BGC.pred))
  for(i in zones){ zone.pred[grep(i,BGC.pred)] <- i }
  assign(paste("zone.change", hist.year, sep="."), zone.pred!=zone) #did the BGC zone change (mapped baseline)
  assign(paste("zone.change.pred", hist.year, sep="."), zone.pred!=zone.pred.ref) #did the BGC zone change (predicted baseline)
  
  print(hist.year)
}

# Future BGC
PredSum <- data.frame()
for(rcp in rcps){
  for(proj.year in proj.years){
    for(GCM in GCMs){
      BGC.pred <- as.character(read.csv(paste("outputs\\BGC.pred",grid, GCM, rcp, proj.year, model,"csv", sep="."), header=F)[,1])
      assign(paste("BGC.pred", GCM, rcp, proj.year, sep="."), BGC.pred) #bgc projection
      PredSum <- rbind(PredSum, as.data.frame(table(BGC.pred)))
      
      #did the BGC change?
      assign(paste("BGC.change", GCM, rcp, proj.year, sep="."), BGC.pred!=BGC) #did the BGC unit change (mapped baseline)
      assign(paste("BGC.change.pred", GCM, rcp, proj.year, sep="."), BGC.pred!=BGC.pred.ref) #did the BGC unit change (predicted baseline)
      
      #did the zone change? 
      zone.pred <- rep(NA, length(BGC.pred))
      for(i in zones){ zone.pred[grep(i,BGC.pred)] <- i }
      assign(paste("zone.change", GCM, rcp, proj.year, sep="."), zone.pred!=zone) #did the BGC zone change (mapped baseline)
      assign(paste("zone.change.pred", GCM, rcp, proj.year, sep="."), zone.pred!=zone.pred.ref) #did the BGC zone change (predicted baseline)
      
      # print(GCM)
    }
    print(proj.year)
  }
  print(rcp)
}

# determine vote winner BGC and ensemble agreement (WARNING: takes about 2 minutes per rcp/proj.year, so i just did RCP45, 2055)
rcp=rcp.focal
proj.year=proj.year.focal
for(proj.year in proj.years){
  temp <- as.data.frame(matrix(rep(NA, length(BGC.pred)*length(GCMs)), nrow=length(BGC.pred), ncol=length(GCMs)))
  for(rcp in rcps){
    for(GCM in GCMs){
      BGC.pred <- get(paste("BGC.pred", GCM, rcp, proj.year, sep="."))
      #add votes to votes matrix
      temp[,which(GCMs==GCM)+length(GCMs)*(which(rcps==rcp)-1)] <- BGC.pred
      print(GCM)
    }
    print(rcp)
  }
  vote.winner <- function(x){return(names(which(table(x)==max(table(x))))[1])}
  agreement <- function(x){return(max(table(x)))}
  assign(paste("BGC.pred.ensemble", proj.year, sep="."), apply(temp, 1, vote.winner))
  # assign(paste("BGC.pred.agreement", rcp, proj.year, sep="."), apply(temp, 1, agreement))
  print(proj.year)
}

# calculate percentage of ensemble projecting an exotic unit
    rcp=rcp.focal
    # for(rcp in rcps){
    proj.year=proj.year.focal
    # for(proj.year in proj.years){
    exotic <- rep(0, length(BGC.pred.ref))
    for(GCM in GCMs){
      BGC.pred <- get(paste("BGC.pred", GCM, rcp, proj.year, sep="."))
      temp <- rep(0, length(BGC.pred.ref))
      temp[-which(BGC.pred%in%unique(BGC))] <- 1
      exotic <- apply(cbind(exotic, temp), 1, sum)
      print(GCM)
    }
    assign(paste("BGC.pred.exotic", rcp, proj.year, sep="."), exotic/length(GCMs))
#     print(proj.year)
#   }
#   print(rcp)
# }



#=============================
## plot of displacement vs MAT change

png(filename=paste("results\\CCISS.BGCEDA.ConsensusEffect.png",sep="."), type="cairo", units="in", width=3.5, height=3.5, pointsize=9, res=300)
par(mar=c(3.25,3.25,0.1,0.1), mgp=c(1.75,0.25,0))

int.y <- 0.07
start.y <- 0.03
pos.y <- seq(start.y,start.y+int.y*3, int.y)
offset.y <- 6
ColScheme=c("dodgerblue", "red")
  proj.year=2055
  # for(proj.year in rev(proj.years)){
    temp.bgc <- rep(NA, length(GCMs))
    temp.zone <- rep(NA, length(GCMs))
    for(rcp in rev(rcps)){
      for(GCM in GCMs){
      temp.bgc[which(GCMs==GCM)+length(GCMs)*(which(rcps==rcp)-1)] <- mean(get(paste("BGC.change.pred", GCM, rcp, proj.year, sep=".")) )
      temp.zone[which(GCMs==GCM)+length(GCMs)*(which(rcps==rcp)-1)] <- mean(get(paste("zone.change.pred", GCM, rcp, proj.year, sep=".")) )
      }
    }
    
    x <- c(get(paste("MAT.change", rcps[1], proj.year, sep=".")), get(paste("MAT.change", rcps[2], proj.year, sep=".")))
    y1 <- temp.bgc
    y2 <- temp.zone
    
    plot(0, xlim=c(0,max(x)+0.5), ylim=c(0,1.01), xaxs="i", yaxs="i", col="white", xaxt="n", yaxt="n",
         xlab=list(bquote(Projected~temperature~change~relative~to~"1961-90"~"("*degree*C*")"), cex=0.8), 
         ylab=list("Projected displacement of modeled BGC unit", cex=0.8))
    axis(1, at=0:8, labels = 0:8, tck=0, cex.axis=0.8)
    axis(2, at=seq(0,1,0.2), labels = paste(seq(0,1,0.2)*100, "%", sep=""), las=2, tck=0, cex.axis=0.8)

        points(x, y1, pch=16)
    points(x, y2, pch=16, col="gray")
    
## add in points for vote winner
BGC.pred.ensemble <- get(paste("BGC.pred.ensemble", proj.year, sep=".")) 
zone.pred.ensemble <- rep(NA, length(BGC.pred.ensemble))
for(i in zones){ zone.pred.ensemble[grep(i,BGC.pred.ensemble)] <- i }

BGC.change.pred.ens <- BGC.pred.ensemble!=BGC.pred.ref
zone.change.pred.ens <- zone.pred.ensemble!=zone.pred.ref

x <- mean(x)
y1 <- mean(BGC.change.pred.ens)
y2 <- mean(zone.change.pred.ens)

points(x, y1, pch=17, cex=2)
points(x, y2, pch=17, cex=2, col="gray")

    # print(proj.year)
    # }
    
legend("bottomright", cex=0.8, legend=c("BGC subzone/variant displacement", "BGC zone displacement", "'Consensus' projection (winning vote)", ""), y.intersp = 1, pch=c(16,16,17, 1), col=c(1,"gray",1,"white"), pt.cex=c(1.2,1.2,1.5, 1), bty="n")
legend("right", cex=0.8, legend=c(proj.year.name[which(proj.years==proj.year)], paste(rcp.name[1], "&", rcp.name[2])), bty="n")

dev.off()
