
##======================================================================================
## CCISS Publication Scripts
## Step 4d - Figures - Suitability for individual species and groups of species
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

#===============================================================================
# Set analysis Parameters
#===============================================================================

setwd("C:\\Colin\\Projects\\2019_CCISS")

grid <- "BC2kmGrid"

GCMs <-  c("ACCESS1-0","CanESM2","CCSM4","CESM1-CAM5","CNRM-CM5","CSIRO-Mk3-6-0", "GFDL-CM3","GISS-E2R", "HadGEM2-ES", "INM-CM4", "IPSL-CM5A-MR", "MIROC-ESM", "MIROC5", "MPI-ESM-LR","MRI-CGCM3")    
rcps <- c("rcp45", "rcp85")
proj.years <- c(2025, 2055, 2085)
edatopes<- c("B2", "C4", "D6")

# Knowledge Tables
treesuit="TreeSpp_ESuit_v11_18"
SiteSeries_Use <-read.csv(paste("InputData/","SiteSeries_Use_5",".csv",sep=""),stringsAsFactors=FALSE,na.strings=".")
spps.lookup <- read.csv("InputData\\Tree speciesand codes_2.0_2May2019.csv")

#===============================================================================
# create a dem from the climateBC input data
#===============================================================================

## create a dem from the climateBC input data
points <- read.csv(paste("InputData\\",grid,".csv", sep=""))
dim(points)

# points <- points[order(points$lon, points$lat),]
# head(points)

# transform latlon points to albers (this only works because the points were originally in albers)
P4S.latlon <- CRS("+proj=longlat +datum=WGS84")
P4S.albers <- CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs")
pts <- points[,4:3]
coordinates(pts) <- pts
projection(pts) <- P4S.latlon
pts.albers <- spTransform(pts, P4S.albers) #reproject the gridpoints

## assemble a simple spatial data frame of the points (seems awkward, but this is how i got the raster to work properly)
pts.df <- as.data.frame(pts.albers)
pts.df<- round(pts.df[,3:4],1)
df = data.frame(z = 1:dim(points)[1],
                xc = pts.df[,1],
                yc = pts.df[,2])
coordinates(df) = ~xc+yc
gridded(df) = TRUE
df = as(df, "SpatialGridDataFrame") # to full grid
plotOrder <- df$z #plotting order relative to input data

## rasterize
dem <- raster(df)
projection(dem) <- P4S.albers
par(mfrow=c(1,1))
plot(dem)
sum(!is.na(values(dem)))

X <- dem
values(X) <- points$el[df$z]
plot(X)

#===============================================================================
# generic spatial data
#===============================================================================
### admin boundaries
bdy.bc <- readOGR("InputData\\BC_AB_US_Shp\\ProvincialOutline.shp")


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

## non-THLB BGCs for exclusion from results
BGC <- points$ID2
BGC <- gsub(" ","",BGC)
BGCs_notin_THLB <- read.csv("InputData\\BGCs_notin_THLB.csv")
BGCs_notin_THLB <- BGCs_notin_THLB$BGC[which(BGCs_notin_THLB$Exlude=="x")]
exclude <- which(BGC%in%BGCs_notin_THLB)

#BGC zones
BGCcolors <- read.csv("C:\\Colin\\Projects\\2019_CCISS\\InputData\\BGCzone_Colorscheme.csv")
zone <- rep(NA, length(BGC))
for(i in BGCcolors$zone){ zone[grep(i,BGC)] <- i }
table(zone)
zone <- factor(zone, levels=c("CDF", "CWH", "MH", "ESSF", "MS", "IDF", "PP", "BG", "ICH", "SBPS", "SBS", "BWBS", "SWB", "CMA", "IMA", "BAFA"))


#===============================================================================
# Manuscript plot for a subset of species
#===============================================================================
edatope.name <- c("Poor-Subxeric", "Medium-Mesic", "Rich-Hygric")
proj.year.name=c("2020s", "2050s", "2080s")
rcp=rcps[1]
proj.year=proj.years[2]

edatope="C4"
# for(edatope in edatopes){

spps.matrix <- matrix(c("Pl", "Fd", "Cw","Ba", "Bl", "Bg", "Yc", "Pa", "Hm", "Lw", "Hw", "Py", "Dr", "Ep", "At"), 5, byrow=T)  

# x11(width=6.5, height=8.5, pointsize=12)

for(matrow in 1:dim(spps.matrix)){
 spps <- spps.matrix[matrow,]

png(filename=paste("Results\\Manu_Suitability_Groups\\CCISS.manu.Suitability",spps[1],spps[2],spps[3], edatope, rcp, proj.year,"png",sep="."), type="cairo", units="in", width=6.5, height=8.5, pointsize=12, res=600)
par(mar=c(0,0,0,0), mfrow=c(3,1), bg="white")

for(spp in spps){

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
  
  ##=================================
  ## initiate plot
  par(plt=c(0,1,0,1), bg=NA)
  plot(0, col="white", xaxt="n", yaxt="n", xlab="", ylab="")
  Common <- as.character(spps.lookup$EnglishName[which(spps.lookup$TreeCode==spp)])
  Latin <- as.character(spps.lookup$ScientificName[which(spps.lookup$TreeCode==spp)])
  panel <- paste("(", LETTERS[which(spps==spp)],")", sep="")
  mtext(if(spp%in%spps.lookup$TreeCode) bquote(bold(.(spp))~"-"~.(Common)) else bquote(bold(.(spp))),
        side=3, line=-1.75, adj=0.01, cex=0.8, font=2)
  # mtext(if(spp%in%spps.lookup$TreeCode) bquote(.(panel)~bold(.(spp))~"-"~.(Common)~"("*italic(.(Latin)*")")) else bquote(.(panel)~bold(.(spp))),
  #       side=3, line=-1.75, adj=0.01, cex=0.8, font=2)
  mtext(paste("Site type: ", edatope, " (", edatope.name[which(edatopes==edatope)], ")", sep=""), side=3, line=-2.75, adj=0.01, cex=0.7, font=1)
  box()
  
  ##=================================
  # map of binary appearance/disappearance
  Suit.ensemble <- as.matrix(ProjSuit)
  Suit.ensemble[Suit.ensemble==5] <- 4
  binary <- rep(0, length(RefSuit))
  binary[outRange.base==T] <- NA
  binary[outRange.base] <- apply(Suit.ensemble[outRange.base,], 1, function(x){return(if((sum(x<4, na.rm=T)/sum(!is.na(x)))>0) sum(x<4, na.rm=T)/sum(!is.na(x)) else NA)})
  binary[outRange.base==F] <- apply(Suit.ensemble[outRange.base==F,], 1, function(x){return(0-sum(x==4, na.rm=T)/sum(!is.na(x)))})
  values(X) <- binary[plotOrder]
  
  breakpoints <- seq(-1,1,0.2); length(breakpoints)
  labels <- c("Retreat", "Expansion")
  ColScheme <- c(brewer.pal(11,"RdBu")[c(1:4)], "grey90", brewer.pal(11,"RdBu")[c(7:11)]); length(ColScheme)
  
  par(plt = c(0.075,0.975,0,1), new = TRUE)
  plot(bdy.bc, border="black", lwd=0.4)
  image(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks=breakpoints, maxpixels= ncell(X))
  # mtext(paste("Edatope:", edatope), side=1, line=-1.5, adj=0.02, cex=1.1, font=2)
  # if(spp==spps[1]){
  xl <- 325000; yb <- 900000; xr <- 400000; yt <- 1525000
  rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme)
  text(rep(xr+10000,length(labels)),seq(yb,yt,(yt-yb)/(length(GCMs)-1))[c(3,9)],labels,pos=4,cex=0.9,font=0.8, srt=90)
  text(rep(xr-20000,length(labels)),seq(yb,yt,(yt-yb)/(length(GCMs)-1))[c(1,8,15)],c("100%", "0%", "100%"),pos=4,cex=0.8,font=1)
  text(xl-30000, mean(c(yb,yt))-30000, paste("Change in presence/absence\n(2050s), % of GCMs"), srt=90, pos=3, cex=0.9, font=2)
  mtext(paste("(", LETTERS[c(2,6,10)][which(spps==spp)],")", sep=""), side=3, line=-2.5, adj=0.22, cex=0.8, font=2)
  # legend("bottomleft", legend=c(spp, paste("Edatope:", edatope), proj.year, rcp, " "), cex=1.4, bty="n", inset=-0.05)
  # }
  # box()
  
  ##=================================
  # historical suitability
  breakseq <- c(0.5,1.5,2.5,3.5,5)
  ColScheme <- c(brewer.pal(9,"Greys")[9], brewer.pal(9,"Greens")[7], brewer.pal(9,"Greens")[4], "white")
  ColScheme <- c("darkgreen", "dodgerblue1", "gold2", "white")
  length(ColScheme)
  
  par(plt = c(0, 0.375, 0, 0.75), new = TRUE)
  values(X) <- RefSuit[plotOrder]
  plot(bdy.bc, border="black", lwd=0.4)
  image(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks=breakseq, maxpixels= ncell(X))
  plot(bdy.bc, add=T, border="black", lwd=0.4)
  # if(spp==spps[1]){
    legend("bottomleft", legend=c("1 (primary)", "2 (secondary)", "3 (tertiary)"), 
         fill=ColScheme, bty="n", cex=0.9, title="Historical suitability", inset=0.015)
  # }
    mtext(paste("(", LETTERS[c(1,5,9)][which(spps==spp)],")", sep=""), side=3, line=-3.75, adj=0.05, cex=0.8, font=2)
    
    ##=================================
    # map of suitability change
    breakpoints <- seq(-3,3,0.5); length(breakpoints)
    labels <- c("-3","-2", "-1", "no change", "+1","+2","+3")
    ColScheme <- c(brewer.pal(11,"RdBu")[c(1,2,3,4,4)], "grey80", brewer.pal(11,"RdBu")[c(7,8,8,9,10,11)]); length(ColScheme)
    
    par(plt = c(0.6, 0.95, 0.3, 1), new = TRUE)
    values(X) <- ChangeSuit.mean[plotOrder]
    plot(bdy.bc, border="black", lwd=0.4)
    image(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks=breakpoints, maxpixels= ncell(X))
    plot(bdy.bc, add=T, border="black", lwd=0.4)
    # if(spp==spps[1]){
    par(xpd=T)
    xl <- 1600000; yb <- 1000000; xr <- 1700000; yt <- 1700000
    rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme)
    text(rep(xr-10000,length(labels)),seq(yb,yt,(yt-yb)/(length(labels)-1)),labels,pos=4,cex=0.8,font=1)
    text(xl-30000, mean(c(yb,yt))-30000, paste("Mean change\nin suitability (", proj.year.name[which(proj.years==proj.year)], ")", sep=""), srt=90, pos=3, cex=0.9, font=2)
    # }
    par(xpd=F)
    mtext(paste("(", LETTERS[c(3,7,11)][which(spps==spp)],")", sep=""), side=3, line=-3.25, adj=0.1, cex=0.8, font=2)
    
    ##=================================
    ## Summary by zone

    # par(mar=c(0,0,0,0), plt = c(0.77, 0.995, 0.001, 0.31), new = TRUE, mgp=c(1.25,0.15,0))
    # plot(0, xlim=c(0,1), ylim=c(0,1), col="white", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
    ChangeSuit.mean[is.na(ChangeSuit.mean)] <- 0
    par(mar=c(4.5,2,0.1,0.1), plt = c(0.79, 0.995, 0.1, 0.3), new = TRUE, mgp=c(1.25,0.15,0))
    ylim=c(-3,3)
    xlim=c(1, length(levels(droplevels(zone))))
    z <- boxplot(ChangeSuit.mean~zone, ylab="", vertical = TRUE, plot=F)
    for(i in 1:length(levels(zone))){ 
      temp <- ChangeSuit.mean[which(zone==levels(zone)[i])]
      z$stats[c(1,5), i] <- quantile(temp[!is.na(temp)],c(0.05, 0.95))
    }
    bxp(z, xlim=xlim, ylim=ylim, xaxt="n", yaxt="n", xaxs="i", ylab="", pch=0,outline=FALSE)
    lines(c(-99,99), c(0,0), lwd=2, col="darkgrey")
    bxp(z, add=T, boxfill = as.character(BGCcolors$HEX[match(levels(zone), BGCcolors$zone)]), xaxt="n", yaxt="n", xaxs="i", ylab="", pch=0,outline=FALSE)
    axis(1, at=1:length(levels(zone)), levels(zone), tick=F, las=2, cex.axis=0.8)
    axis(2,at=seq(ylim[1], ylim[2], 3), seq(ylim[1], ylim[2], 3), las=2, tck=0)
    mtext(paste("(", LETTERS[c(4,8,12)][which(spps==spp)],")", sep=""), side=3, line=1, adj=0.975, cex=0.8, font=2)
    mtext("Mean change in suitability", side=3, line=0.1, adj=.975, cex=0.55, font=2)

    print(spp)
}
dev.off()
print(matrow)
}

#===============================================================================
# Manuscript plot of all edatopes for one species
#===============================================================================
rcp=rcps[1]
proj.year=proj.years[2]

# select the species to run the analysis on
spps <- unique(S1$Spp)
spps <- spps[-which(spps=="X")]
spps.candidate <- spps.lookup$TreeCode[-which(spps.lookup$Exclude=="x")]
spps <- spps[which(spps%in%spps.candidate)] 

spp="Sx"

for(spp in spps){

  png(filename=paste("Results\\Manu_Suitability_Indiv\\CCISS.manu.Suitability.edatopes",spp, rcp, proj.year,"png",sep="."), type="cairo", units="in", width=6.5, height=8.5, pointsize=12, res=600)
  par(mar=c(0,0,0,0), mfrow=c(3,1), bg="white")
  
  for(edatope in edatopes){
    
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
    
    ##=================================
    ## initiate plot
    par(plt=c(0,1,0,1), bg=NA)
    plot(0, col="white", xaxt="n", yaxt="n", xlab="", ylab="")
    Common <- as.character(spps.lookup$EnglishName[which(spps.lookup$TreeCode==spp)])
    Latin <- as.character(spps.lookup$ScientificName[which(spps.lookup$TreeCode==spp)])
    panel <- paste("(", LETTERS[which(spps==spp)],")", sep="")
    mtext(paste("Site type: ", edatope, " (", edatope.name[which(edatopes==edatope)], ")", sep=""), side=3, line=-1.55, adj=0.01, cex=0.8, font=2)
    mtext(if(spp%in%spps.lookup$TreeCode) bquote(bold(.(spp))~"-"~.(Common)) else bquote(bold(.(spp))), side=3, line=-2.75, adj=0.01, cex=0.7, font=1)
    box()
    
    ##=================================
    # map of binary appearance/disappearance
    Suit.ensemble <- as.matrix(ProjSuit)
    Suit.ensemble[Suit.ensemble==5] <- 4
    binary <- rep(0, length(RefSuit))
    binary[outRange.base==T] <- NA
    binary[outRange.base] <- apply(Suit.ensemble[outRange.base,], 1, function(x){return(if((sum(x<4, na.rm=T)/sum(!is.na(x)))>0) sum(x<4, na.rm=T)/sum(!is.na(x)) else NA)})
    binary[outRange.base==F] <- apply(Suit.ensemble[outRange.base==F,], 1, function(x){return(0-sum(x==4, na.rm=T)/sum(!is.na(x)))})
    values(X) <- binary[plotOrder]
    
    breakpoints <- seq(-1,1,0.2); length(breakpoints)
    labels <- c("Retreat", "Expansion")
    ColScheme <- c(brewer.pal(11,"RdBu")[c(1:4)], "grey90", brewer.pal(11,"RdBu")[c(7:11)]); length(ColScheme)
    
    par(plt = c(0.075,0.975,0,1), new = TRUE)
    plot(bdy.bc, border="black", lwd=0.4)
    image(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks=breakpoints, maxpixels= ncell(X))
    # mtext(paste("Edatope:", edatope), side=1, line=-1.5, adj=0.02, cex=1.1, font=2)
    # if(spp==spps[1]){
    xl <- 325000; yb <- 900000; xr <- 400000; yt <- 1525000
    rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme)
    text(rep(xr+10000,length(labels)),seq(yb,yt,(yt-yb)/(length(GCMs)-1))[c(3,9)],labels,pos=4,cex=0.9,font=0.8, srt=90)
    text(rep(xr-20000,length(labels)),seq(yb,yt,(yt-yb)/(length(GCMs)-1))[c(1,8,15)],c("100%", "0%", "100%"),pos=4,cex=0.8,font=1)
    text(xl-30000, mean(c(yb,yt))-30000, paste("Change in presence/absence\n(2050s), % of GCMs"), srt=90, pos=3, cex=0.9, font=2)
    mtext(paste("(", LETTERS[c(2,6,10)][which(edatopes==edatope)],")", sep=""), side=3, line=-2.75, adj=0.22, cex=0.8, font=2)
    # legend("bottomleft", legend=c(spp, paste("Edatope:", edatope), proj.year, rcp, " "), cex=1.4, bty="n", inset=-0.05)
    # }
    # box()
    
    ##=================================
    # historical suitability
    breakseq <- c(0.5,1.5,2.5,3.5,5)
    ColScheme <- c(brewer.pal(9,"Greys")[9], brewer.pal(9,"Greens")[7], brewer.pal(9,"Greens")[4], "white")
    ColScheme <- c("darkgreen", "dodgerblue1", "gold2", "white")
    length(ColScheme)
    
    par(plt = c(0, 0.375, 0, 0.75), new = TRUE)
    values(X) <- RefSuit[plotOrder]
    plot(bdy.bc, border="black", lwd=0.4)
    image(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks=breakseq, maxpixels= ncell(X))
    plot(bdy.bc, add=T, border="black", lwd=0.4)
    # if(spp==spps[1]){
    legend("bottomleft", legend=c("1 (primary)", "2 (secondary)", "3 (tertiary)"), 
           fill=ColScheme, bty="n", cex=0.9, title="Historical suitability", inset=0.015)
    # }
    mtext(paste("(", LETTERS[c(1,5,9)][which(edatopes==edatope)],")", sep=""), side=3, line=-3.75, adj=0.05, cex=0.8, font=2)
    
    ##=================================
    # map of suitability change
    breakpoints <- seq(-3,3,0.5); length(breakpoints)
    labels <- c("-3","-2", "-1", "no change", "+1","+2","+3")
    ColScheme <- c(brewer.pal(11,"RdBu")[c(1,2,3,4,4)], "grey80", brewer.pal(11,"RdBu")[c(7,8,8,9,10,11)]); length(ColScheme)
    
    par(plt = c(0.6, 0.95, 0.3, 1), new = TRUE)
    values(X) <- ChangeSuit.mean[plotOrder]
    plot(bdy.bc, border="black", lwd=0.4)
    image(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks=breakpoints, maxpixels= ncell(X))
    plot(bdy.bc, add=T, border="black", lwd=0.4)
    # if(spp==spps[1]){
    par(xpd=T)
    xl <- 1600000; yb <- 1000000; xr <- 1700000; yt <- 1700000
    rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme)
    text(rep(xr-10000,length(labels)),seq(yb,yt,(yt-yb)/(length(labels)-1)),labels,pos=4,cex=0.8,font=1)
    text(xl-30000, mean(c(yb,yt))-30000, paste("Mean change\nin suitability (", proj.year.name[which(proj.years==proj.year)], ")", sep=""), srt=90, pos=3, cex=0.9, font=2)
    # }
    par(xpd=F)
    mtext(paste("(", LETTERS[c(3,7,11)][which(edatopes==edatope)],")", sep=""), side=3, line=-3.25, adj=0.1, cex=0.8, font=2)
    
    ##=================================
    ## Summary by zone
    
    # par(mar=c(0,0,0,0), plt = c(0.77, 0.995, 0.001, 0.31), new = TRUE, mgp=c(1.25,0.15,0))
    # plot(0, xlim=c(0,1), ylim=c(0,1), col="white", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
    ChangeSuit.mean[is.na(ChangeSuit.mean)] <- 0
    par(mar=c(4.5,2,0.1,0.1), plt = c(0.79, 0.995, 0.1, 0.3), new = TRUE, mgp=c(1.25,0.15,0))
    ylim=c(-3,3)
    xlim=c(1, length(levels(droplevels(zone))))
    z <- boxplot(ChangeSuit.mean~zone, ylab="", vertical = TRUE, plot=F)
    for(i in 1:length(levels(zone))){ 
      temp <- ChangeSuit.mean[which(zone==levels(zone)[i])]
      z$stats[c(1,5), i] <- quantile(temp[!is.na(temp)],c(0.05, 0.95))
    }
    bxp(z, xlim=xlim, ylim=ylim, xaxt="n", yaxt="n", xaxs="i", ylab="", pch=0,outline=FALSE)
    lines(c(-99,99), c(0,0), lwd=2, col="darkgrey")
    bxp(z, add=T, boxfill = as.character(BGCcolors$HEX[match(levels(zone), BGCcolors$zone)]), xaxt="n", yaxt="n", xaxs="i", ylab="", pch=0,outline=FALSE)
    axis(1, at=1:length(levels(zone)), levels(zone), tick=F, las=2, cex.axis=0.8)
    axis(2,at=seq(ylim[1], ylim[2], 3), seq(ylim[1], ylim[2], 3), las=2, tck=0)
    mtext(paste("(", LETTERS[c(4,8,12)][which(edatopes==edatope)],")", sep=""), side=3, line=1, adj=0.975, cex=0.8, font=2)
    mtext("Mean change in suitability", side=3, line=0.1, adj=.975, cex=0.55, font=2)
    
    print(edatope)
  }
  dev.off()
  print(spp)
}

  
  
#===============================================================================
# Basic plot for a single species
#===============================================================================
# spps.name <- c("Yellow cedar", "Interior spruce", "Ponderosa pine", "Lodgepole pine", "Western larch", "Western hemlock", "Mountain hemlock", "Douglas-fir", "Red alder", "Western redcedar", "Subalpine fir", "Grand fir")
edatope.name <- c("Subxeric-poor", "Mesic-medium", "Hygric-rich")
proj.year.name=c("2020s", "2050s", "2080s")
rcp=rcps[1]
proj.year=proj.years[1]

edatope="C4"
for(proj.year in proj.years){
  for(edatope in edatopes){
    
    for(spp in spps){
      
      png(filename=paste("Results\\SI_Suitability_Basic\\CCISS.manu.Suitability.basic",spp, edatope, rcp, proj.year,"png",sep="."), type="cairo", units="in", width=6.5, height=6, pointsize=10, res=600)
      par(mar=c(0,0,0,0), mfrow=c(2,2), bg="white")
      
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
      
      values(X) <- RefSuit[plotOrder]
      plot(bdy.bc, border="black", lwd=0.4)
      image(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks=breakseq, maxpixels= ncell(X))
      plot(bdy.bc, add=T, border="black", lwd=0.4)
      # if(spp==spps[1]){
      legend("topright", legend=c("1 (primary)", "2 (secondary)", "3 (tertiary)"),
             fill=ColScheme, bty="n", cex=0.9, title="Historical suitability", inset=0.015)
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
      
      # map of suitability change
      breakpoints <- seq(-3,3,0.5); length(breakpoints)
      labels <- c("-3","-2", "-1", "no change", "+1","+2","+3")
      ColScheme <- c(brewer.pal(11,"RdBu")[c(1,2,3,4,4)], "grey80", brewer.pal(11,"RdBu")[c(7,8,8,9,10,11)]); length(ColScheme)
      
      values(X) <- ChangeSuit.mean[plotOrder]
      plot(bdy.bc, border="black", lwd=0.4)
      image(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks=breakpoints, maxpixels= ncell(X))
      plot(bdy.bc, add=T, border="black", lwd=0.4)
      # if(spp==spps[1]){
      par(xpd=T)
      xl <- 1600000; yb <- 1000000; xr <- 1700000; yt <- 1700000
      rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme)
      text(rep(xr-10000,length(labels)),seq(yb,yt,(yt-yb)/(length(labels)-1)),labels,pos=4,cex=0.8,font=1)
      text(xl-30000, mean(c(yb,yt))-30000, paste("Mean change\nin suitability (", proj.year.name[which(proj.years==proj.year)], ")", sep=""), srt=90, pos=3, cex=0.9, font=2)
      # }
      par(xpd=F)
      
      
      # map of ensemble agreement on trend
      increasing <- which(ChangeSuit.mean>0)
      decreasing <- which(ChangeSuit.mean<0)
      bifurc <- rep(NA, length(ChangeSuit.mean))
      bifurc[outRange==F] <- 0
      bifurc[increasing] <- apply(ChangeSuit[increasing,], 1, function(x){return(sum(x< 0, na.rm=T)/length(x))})
      bifurc[decreasing] <- apply(ChangeSuit[decreasing,], 1, function(x){return(sum(x> 0, na.rm=T)/length(x))})
      values(X) <- bifurc[plotOrder]
      
      breakpoints <- seq(0,0.5,0.1); length(breakpoints)
      labels <- c("High", "Medium", "Low")
      # ColScheme <- c(brewer.pal(11,"RdBu")[c(4,2,1)], "black", brewer.pal(11,"RdBu")[c(8,10,11)]); length(ColScheme)
      # ColScheme <- brewer.pal(11,"RdBu")[-c(2,5,7,10)]; length(ColScheme)
      ColScheme <- c("grey80", brewer.pal(9,"YlOrRd")[c(3,5,7,9)]); length(ColScheme)
      
      plot(bdy.bc, border="black", lwd=0.4)
      image(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks=breakpoints, maxpixels= ncell(X))
      # mtext(paste("Edatope:", edatope), side=1, line=-1.5, adj=0.02, cex=1.1, font=2)
      # if(spp==spps[1]){
      xl <- 1600000; yb <- 1000000; xr <- 1700000; yt <- 1700000
      rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme)
      text(rep(xr,length(labels)),seq(yb,yt,(yt-yb)/(length(labels)*5-1))[c(3,13)],labels[c(1,3)],pos=4,cex=0.9,font=1)
      text(rep(xr,length(labels)),seq(yb,yt,(yt-yb)/(length(labels)*5-1))[c(1,15)],paste(round(breakpoints*100), "%", sep="")[c(1, length(breakpoints))],pos=4,cex=0.9,font=1)
      text(xl-30000, mean(c(yb,yt))-30000, paste("Ensemble agreement\non trend"), srt=90, pos=3, cex=1, font=2)
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
      values(X) <- binary[plotOrder]
      
      breakpoints <- seq(-1,1,0.2); length(breakpoints)
      labels <- c("Retreat", "Expansion")
      ColScheme <- c(brewer.pal(11,"RdBu")[c(1:4)], "grey90", brewer.pal(11,"RdBu")[c(7:11)]); length(ColScheme)
      
      plot(bdy.bc, border="black", lwd=0.4)
      image(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks=breakpoints, maxpixels= ncell(X))
      # mtext(paste("Edatope:", edatope), side=1, line=-1.5, adj=0.02, cex=1.1, font=2)
      # if(spp==spps[1]){
      xl <- 1600000; yb <- 900000; xr <- 1700000; yt <- 1700000
      rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme)
      text(rep(xr+10000,length(labels)),seq(yb,yt,(yt-yb)/(length(GCMs)-1))[c(3,9)],labels,pos=4,cex=0.9,font=0.8, srt=90)
      text(rep(xr-20000,length(labels)),seq(yb,yt,(yt-yb)/(length(GCMs)-1))[c(1,8,15)],c("100%", "0%", "100%"),pos=4,cex=0.8,font=1)
      text(xl-30000, mean(c(yb,yt))-30000, paste("Change in presence/absence\n(% of models)"), srt=90, pos=3, cex=0.9, font=2)
      # legend("bottomleft", legend=c(spp, paste("Edatope:", edatope), proj.year, rcp, " "), cex=1.4, bty="n", inset=-0.05)
      # }
      # box()
      
      print(spp)
      dev.off()
      
    }
    print(edatope)
  }
  print(proj.year)
}

