
##======================================================================================
## CCISS Publication Scripts
## Step 4f - Figures - Manuscript plot of Suitability Richness
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
hist.years <- c(1995, 2005, 2014, 2017)
proj.years <- c(2025, 2055, 2085)
edatopes<- c("B2", "C4", "D6")
edatope.name <- c("nutrient poor, subxeric", "nutrient-medium, mesic", "nutrient-rich, hygric")

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
# generate the vector of predicted referencc BGCs 
#===============================================================================

BGC <- as.character(read.csv(paste("OutputData\\BGC.pred", grid, "ref.csv", sep="."))[,1])

#BGC zones
BGCcolors <- read.csv("C:\\Colin\\Projects\\2019_CCISS\\InputData\\BGCzone_Colorscheme.csv")
zone <- rep(NA, length(BGC))
for(i in BGCcolors$zone){ zone[grep(i,BGC)] <- i }
table(zone)
zone <- factor(zone, levels=c("CDF", "CWH", "MH", "ESSF", "MS", "IDF", "PP", "BG", "ICH", "SBPS", "SBS", "BWBS", "SWB", "CMA", "IMA", "BAFA"))


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

#===============================================================================
# calculate mean optionality and turnover for each edatope/year/rcp
#===============================================================================

# read in metrics for reference period and historical decades
for(edatope in edatopes){
  assign(paste("SuitRichness.ref", edatope, sep="."), read.csv(paste("OutputData\\SuitRichness.ref", grid, edatope, "csv", sep="."))[,1])
  assign(paste("SppRichness.ref", edatope, sep="."), read.csv(paste("OutputData\\SppRichness.ref", grid, edatope, "csv", sep="."))[,1])
  print(edatope)
}

for(edatope in edatopes){
  for(rcp in rcps){
    for(proj.year in proj.years){
      assign(paste("SuitRichness", rcp, proj.year, edatope, sep="."), read.csv(paste("OutputData\\SuitRichness", grid, rcp, proj.year, edatope, "csv", sep="."))[,1])
      assign(paste("SuitRichnessChange", rcp, proj.year, edatope, sep="."), read.csv(paste("OutputData\\SuitRichnessChange", grid, rcp, proj.year, edatope, "csv", sep="."))[,1])
      assign(paste("SppRichness", rcp, proj.year, edatope, sep="."), read.csv(paste("OutputData\\SppRichness", grid, rcp, proj.year, edatope, "csv", sep="."))[,1])
      assign(paste("SppRichnessChange", rcp, proj.year, edatope, sep="."), read.csv(paste("OutputData\\SppRichnessChange", grid, rcp, proj.year, edatope, "csv", sep="."))[,1])
      assign(paste("SuitTurnover", rcp, proj.year, edatope, sep="."), read.csv(paste("OutputData\\SuitTurnover", grid, rcp, proj.year, edatope, "csv", sep="."))[,1])
      assign(paste("SppTurnover", rcp, proj.year, edatope, sep="."), read.csv(paste("OutputData\\SppTurnover", grid, rcp, proj.year, edatope, "csv", sep="."))[,1])
      assign(paste("SuitPersistence", rcp, proj.year, edatope, sep="."), read.csv(paste("OutputData\\SuitPersistence", grid, rcp, proj.year, edatope, "csv", sep="."))[,1])
      assign(paste("SppPersistence", rcp, proj.year, edatope, sep="."), read.csv(paste("OutputData\\SppPersistence", grid, rcp, proj.year, edatope, "csv", sep="."))[,1])
      print(proj.year)
    }
    print(rcp)
  }
  print(edatope)
}



#===============================================================================
# BGC projections
#===============================================================================

# Future BGC
PredSum <- data.frame()
for(rcp in rcps){
  for(proj.year in proj.years){
    for(GCM in GCMs){
      BGC.pred <- as.character(read.csv(paste("OutputData\\BGC.pred",grid, GCM, rcp, proj.year,".csv", sep=""))[,1])
      assign(paste("BGC.pred", GCM, rcp, proj.year, sep="."), BGC.pred) #bgc projection
      PredSum <- rbind(PredSum, as.data.frame(table(BGC.pred)))
      # print(GCM)
    }
    print(proj.year)
  }
  print(rcp)
}




 ################################
 ## Manuscript figure
 ################################
edatope=edatopes[2]
rcp=rcps[1]
proj.year=proj.years[1]
for(proj.year in proj.years){
  for(edatope in edatopes){

 #===============================================================================
 # assemble the species mix for the reference period
 #===============================================================================
 # for(edatope in edatopes){
 comm.ref <- as.data.frame(matrix(rep(NA, dim(points)[1]*length(spps)), dim(points)[1], length(spps)))
 for(spp in spps){
   Suit <- read.csv(paste("OutputData\\Suit.ref", grid, spp, edatope, "csv", sep="."))[,1]
   Suit[is.na(Suit)] <- 5  #XXX note this is different from the equivalent line for the other time periods. 
   Suit <- 1-(Suit-1)/4
   comm.ref[,which(spps==spp)] <- Suit
 }
 names(comm.ref) <- spps
 assign(paste("comm.ref", edatope, sep="."), comm.ref)
 # Optionality <- apply(comm.ref, 1, sum)
 # assign(paste("Optionality.ref", edatope, sep=""), Optionality)
 # write.csv(Optionality, paste("OutputData\\Optionality.ref", grid, edatope, "csv", sep="."), row.names = F)
 #   print(edatope)
 # }
 
 # ===============================================================================
 # assemble the projected community
 # ===============================================================================
 
  for(GCM in GCMs){
   # for(edatope in edatopes){
   #   for(rcp in rcps){
   #     for(proj.year in proj.years){
   
   comm.proj <- as.data.frame(matrix(rep(NA, dim(points)[1]*length(spps)), dim(points)[1], length(spps)))
   for(spp in spps){
     Suit <- read.csv(paste("OutputData\\Suit", grid, GCM, rcp, proj.year, spp, edatope, "csv", sep="."))[,1]
     Suit[Suit==5] <- 5
     Suit <- 1-(Suit-1)/4
     comm.proj[,which(spps==spp)] <- Suit
   }
   names(comm.proj) <- spps
  assign(paste("comm.proj", GCM, sep="."), comm.proj)
   
   #        # print(proj.year)
   #     }
   #     # print(rcp)
   #   }
   #   print(edatope)
   # }
   # print(GCM)
 }
 
 # ===============================================================================
 # get metrics for the rcp, time period and edatope
 # ===============================================================================
 
 SuitRichness.ref <- get(paste("SuitRichness.ref", edatope, sep="."))
 SuitRichness.proj <- get(paste("SuitRichness", rcp, proj.year, edatope, sep="."))
 SuitRichnessChange <- get(paste("SuitRichnessChange", rcp, proj.year, edatope, sep="."))
 SuitRichnessChangePct <- SuitRichness.proj/SuitRichness.ref
 SuitRichnessChangePct[!is.finite(SuitRichnessChangePct)] <- NA
 SppRichness.ref <- get(paste("SppRichness.ref", edatope, sep="."))
 SppRichness.proj <- get(paste("SppRichness", rcp, proj.year, edatope, sep="."))
 SppRichnessChange <- get(paste("SppRichnessChange", rcp, proj.year, edatope, sep="."))
 SppRichnessChangePct <- SppRichness.proj/SppRichness.ref
 SppRichnessChangePct[!is.finite(SppRichnessChangePct)] <- NA
 SuitTurnover <- get(paste("SuitTurnover", rcp, proj.year, edatope, sep="."))
 SppTurnover <- get(paste("SppTurnover", rcp, proj.year, edatope, sep="."))
 SuitTurnover[which(SuitTurnover==99)] <- NA
 SppTurnover[which(SppTurnover==99)] <- NA
 SuitPersistence <- get(paste("SuitPersistence", rcp, proj.year, edatope, sep="."))
 SppPersistence <- get(paste("SppPersistence", rcp, proj.year, edatope, sep="."))
 
 # ===============================================================================
 # calculate richness with truncation of ensemble mean suitability at 0.375 (treat these as unsuitable)
 # ===============================================================================
 
 comm.proj.ensMean <- as.data.frame(matrix(rep(0, dim(points)[1]*length(spps)), dim(points)[1], length(spps)))
 
 for(GCM in GCMs){
   comm.proj.ensMean <- comm.proj.ensMean+get(paste("comm.proj", GCM, sep="."))
    # print(GCM)
   }
 comm.proj.ensMean <- comm.proj.ensMean/length(GCMs)
 names(comm.proj.ensMean) <- spps
 head(comm.proj.ensMean)

 comm.proj.ensMean.trunc <- as.matrix(comm.proj.ensMean)
 comm.proj.ensMean.trunc[which(comm.proj.ensMean.trunc<0.375)] <- 0
 comm.proj.ensMean.trunc <- as.data.frame(comm.proj.ensMean.trunc)
 head(comm.proj.ensMean.trunc)
 SuitRichness.proj.truncate <- apply(comm.proj.ensMean.trunc, 1, sum)
 SuitRichnessChangePct.truncate <- SuitRichness.proj.truncate/SuitRichness.ref
 SuitRichnessChangePct.truncate[!is.finite(SuitRichnessChangePct.truncate)] <- NA
 

 metric <- "SuitRichnessChangePct"
 
 
 # x11(width=6.5, height=5, pointsize=8)
 png(filename=paste("Results\\Manu_Richness\\CCISS_manu_", metric, edatope, rcp, proj.year,"png",sep="."), type="cairo", units="in", width=6.5, height=5, pointsize=8, res=400)
 # pdf(file=paste("Results\\CCISS_SummaryByBGC_", metric,".pdf",sep=""),  width=7.5, height=5.625, pointsize=15)

 par(mar=c(0,0,0,0))
 plot(0, xlim=c(0,1), ylim=c(0,1), col="white", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
 # box()
 
 #map
 par(mar=c(0,0,0,0), plt = c(0.1, 0.9, 0.125, 0.999), new = TRUE)
 # SuitTurnover <- get(paste("SuitTurnover", rcp, proj.year, edatope, sep="."))
 # SppTurnover <- get(paste("SppTurnover", rcp, proj.year, edatope, sep="."))
 
 ylim <- c(-3,3)
 y <- get(metric)
 y[y<2^(ylim[1])] <- 2^(ylim[1])
 y <- log2(y)
 values(X) <- y[plotOrder]
 
 breakpoints <- c(-99,seq(-2,2,0.5),99); length(breakpoints)
 # breakpoints <- c(-99,-55,seq(-2,2,0.25),55,99); length(breakpoints)
 labels <- c(paste(round(2^(breakpoints[c(2,4)])*100),"%", sep=""), "no change", paste(round(2^(breakpoints[length(breakpoints)-c(3,1)])*100),"%", sep=""))
 ColScheme <- brewer.pal(11,"RdBu")[-6]; length(ColScheme)
 # ColScheme <- c(rep(brewer.pal(11,"RdBu")[1:4], each=2), brewer.pal(11,"RdBu")[5],rep("gray90", 2),brewer.pal(11,"RdBu")[7], rep(brewer.pal(11,"RdBu")[8:11], each=2)); length(ColScheme)
 
 plot(bdy.bc, border="black", lwd=0.4)
 image(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks=breakpoints, maxpixels= ncell(X))
 # mtext(paste("Edatope:", edatope), side=1, line=-1.5, adj=0.02, cex=1.1, font=2)
 par(xpd=T)
 xl <- 150000; yb <- 1100000; xr <- 250000; yt <- 1550000
 xl <- 2025000; yb <- 500000; xr <- 2100000; yt <- 950000
 rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme)
 rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme)
 text(rep(xl,length(labels)),seq(yb,yt,(yt-yb)/(length(breakpoints)-1))[c(2,4,6,8,10)],labels,pos=2,cex=1,font=1)
 # text(rep(xl,2),c(yb,yt)+c(-10000, 10000),c("No change", "Full turnover"),pos=c(1,3),cex=1,font=1)
 text(xr+40000, mean(c(yb,yt))-30000, paste("Relative suitability richness"), srt=90, pos=3, cex=1, font=2)
 # rect(xl,  yt+20000,  xr,  yt+60000,  col=ColScheme[length(ColScheme)])
 # text(xr,  yt+40000,  bquote(">"*.(breakseq[3])*sigma),pos=4,cex=1,font=1)  
 par(xpd=F)
 mtext("(A)", side=3, line=-3.5, adj=0.025, cex=1, font=2)
 mtext(paste(edatope, " edatope (", edatope.name[which(edatopes==edatope)], " sites)", sep=""), side=3, line=-1.5, adj=0.3, cex=1, font=2)
 

 #select a single gridpoint within a specified BGC
 focal.bgcs <- c("ICHxw", "SBSwk1", "ESSFmv4", "SBSmc2")
 q.selects <- c(0.885, 0.5, 0.995, 0.975) #quantile of turnover to select the example grid point from 
 line.x <- c(0, 300000, 300000, -400000)
 line.y <- c(-200000, 300000, 500000, 300000)
 pts.select <- c(106318, 130998, 179068, 53417) # use this if you have already selected the points
 for(i in 1:4){
 focal.bgc <- which(BGC==focal.bgcs[i])
 # pt <- focal.bgc[which(SuitTurnover[focal.bgc]==min(SuitTurnover[focal.bgc][SuitTurnover[focal.bgc]>quantile(SuitTurnover[focal.bgc], q.selects[i], na.rm=T)], na.rm=T))[1]] #select the grid point
 # pt <- pt[1]
 pt <- pts.select[i]
  lines(c(pts.df[pt,1], pts.df[pt,1]+line.x[i]), c(pts.df[pt,2], pts.df[pt,2]+line.y[i]), lwd=1.5)
 points(pts.df[pt,], pch=21, cex=1.75, lwd=1.5, bg=ColScheme[cut(y[pt],breaks=breakpoints)])
 # print(y[pt])
 # print(pt)
 }
 
 ## Callout box for Location 1
 i=1
 bgc.select <- focal.bgcs[i]
 focal.bgc <- which(BGC==focal.bgcs[i])
 pt <- pts.select[i] 

 #convert fractional suitability to standard suitability
 comm.ref.pt <- comm.ref[pt,]
 comm.ref.pt <- round((1-comm.ref.pt)*4+1)
 comm.ref.pt[comm.ref.pt==4] <- 5
 
 #calculate ensemble mean fractional suitability and convert to standard suitability
 for(GCM in GCMs){
   if(GCM==GCMs[1]) comm.proj.pt <- get(paste("comm.proj", GCM, sep="."))[pt,] else comm.proj.pt <- rbind(comm.proj.pt,get(paste("comm.proj", GCM, sep="."))[pt,])
 }
 comm.proj.pt.ensMean <- apply(comm.proj.pt, 2, FUN=mean, na.rm=T)
 comm.proj.pt.ensMean <- round((1-comm.proj.pt.ensMean)*4+1)
 comm.proj.pt.ensMean[comm.proj.pt.ensMean==4] <- 5
 
 std.ref <- paste(paste(names(comm.ref.pt)[which(comm.ref.pt==1)],collapse=""), "(", 
                  paste(names(comm.ref.pt)[which(comm.ref.pt==2)],collapse=""), ")((", 
                  paste(names(comm.ref.pt)[which(comm.ref.pt==3)],collapse=""), "))", sep="")
 std.proj <- paste(paste(names(comm.proj.pt.ensMean)[which(comm.proj.pt.ensMean==1)],collapse=""), "(", 
                   paste(names(comm.proj.pt.ensMean)[which(comm.proj.pt.ensMean==2)],collapse=""), ")((", 
                   paste(names(comm.proj.pt.ensMean)[which(comm.proj.pt.ensMean==3)],collapse=""), "))", sep="")
 
 BGC.list <- vector()
 for(GCM in GCMs){
   BGC.list[which(GCMs==GCM)] <- get(paste("BGC.pred", GCM, rcp, proj.year, sep="."))[pt]
   # print(BGC.pred[pt])
   # print(GCM)
 }
 BGC.list 
 
 par(mar=c(0,0,0,0), plt = c(0.6, 0.99, 0.01, 0.175), new = TRUE, mgp=c(2,0.1,0))
 plot(0, xlim=c(0,1), ylim=c(0,1), col="white", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
 rect(-9,-9,9,9, col="white", lwd=0)
 mtext("(E)", side=3, line=-1, adj=-0.08, cex=1, font=2)
 box()
 
 par(mar=c(0,0,0,0), plt = c(0.6, 0.99, 0.01, 0.175), new = TRUE, mgp=c(2,0.1,0))
 plot(0, xlim=c(0,1), ylim=c(0,1), col="white", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
 text(-0.025, .95, paste("Historical bioclimate"), pos=4, cex=0.8, font=2, offset=0.1)
 text(0, 0.825, bgc.select, pos=4, cex=0.8, offset=0.1)
 text(-.025, 0.65, "Historical suitabilities", pos=4, cex=0.8, font=2, offset=0.1)
 text(0, 0.525, std.ref, pos=4, cex=0.8, offset=0.1)
 text(-.025, 0.35, paste("Mean projected suitabilities", sep=""), pos=4, cex=0.8, font=2, offset=0.1)
 text(0, 0.225, std.proj, pos=4, cex=0.8, offset=0.1)
 text(-.025, 0.05, paste("Mean rel. suit. richness: ", round(get(metric)[pt]*100), "%", sep=""), pos=4, cex=0.8, offset=0.1, font=2)
 text(0.605, .95, paste("Projected bioclimates"), pos=4, cex=0.8, font=2, offset=0.1)

  par(mar=c(3,0,1,0), plt = c(0.85, 0.985, 0.07, 0.145), new = TRUE, mgp=c(0.5,0.1,0))
  bar.bgc <- names(rev(sort(table(BGC.list))))
  bar.zone <- rep(NA, length(bar.bgc))
  for(i in BGCcolors$zone){ bar.zone[grep(i,bar.bgc)] <- i }
  barplot(rev(sort(table(BGC.list))), col=as.character(BGCcolors$HEX[match(bar.zone, BGCcolors$zone)]), horiz=F, las=2, ylab=list("# GCMs", cex=0.7), yaxt="n", cex.names=0.7)
  par(mgp=c(1,0,0))
 axis(2, at=seq(0,12,4), labels=seq(0,12,4), las=1, tck=0, cex.axis=0.8, lwd=0)
 
 
 ## Callout box for Location 2
 i=2
 bgc.select <- focal.bgcs[i]
 focal.bgc <- which(BGC==focal.bgcs[i])
 pt <- pts.select[i] 

 #convert fractional suitability to standard suitability
 comm.ref.pt <- comm.ref[pt,]
 comm.ref.pt <- round((1-comm.ref.pt)*4+1)
 comm.ref.pt[comm.ref.pt==4] <- 5
 
 #calculate ensemble mean fractional suitability and convert to standard suitability
 for(GCM in GCMs){
   if(GCM==GCMs[1]) comm.proj.pt <- get(paste("comm.proj", GCM, sep="."))[pt,] else comm.proj.pt <- rbind(comm.proj.pt,get(paste("comm.proj", GCM, sep="."))[pt,])
 }
 comm.proj.pt.ensMean <- apply(comm.proj.pt, 2, FUN=mean, na.rm=T)
 comm.proj.pt.ensMean <- round((1-comm.proj.pt.ensMean)*4+1)
 comm.proj.pt.ensMean[comm.proj.pt.ensMean==4] <- 5
 
 std.ref <- paste(paste(names(comm.ref.pt)[which(comm.ref.pt==1)],collapse=""), "(", 
                  paste(names(comm.ref.pt)[which(comm.ref.pt==2)],collapse=""), ")((", 
                  paste(names(comm.ref.pt)[which(comm.ref.pt==3)],collapse=""), "))", sep="")
 std.proj <- paste(paste(names(comm.proj.pt.ensMean)[which(comm.proj.pt.ensMean==1)],collapse=""), "(", 
                   paste(names(comm.proj.pt.ensMean)[which(comm.proj.pt.ensMean==2)],collapse=""), ")((", 
                   paste(names(comm.proj.pt.ensMean)[which(comm.proj.pt.ensMean==3)],collapse=""), "))", sep="")
 
 BGC.list <- vector()
 for(GCM in GCMs){
   BGC.list[which(GCMs==GCM)] <- get(paste("BGC.pred", GCM, rcp, proj.year, sep="."))[pt]
   # print(BGC.pred[pt])
   # print(GCM)
 }
 BGC.list 
 
 plt.callout <- c(0.65, 0.99, 0.575, 0.755)
 
 par(mar=c(0,0,0,0), plt = plt.callout, new = TRUE, mgp=c(2,0.1,0))
 plot(0, xlim=c(0,1), ylim=c(0,1), col="white", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
 rect(-9,-9,9,9, col="white", lwd=0)
 mtext("(D)", side=3, line=0.25, adj=0.01, cex=1, font=2)
 box()
 
 par(mar=c(0,0,0,0), plt = plt.callout, new = TRUE, mgp=c(2,0.1,0))
 plot(0, xlim=c(0,1), ylim=c(0,1), col="white", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
 text(-0.025, .95, paste("Historical bioclimate"), pos=4, cex=0.8, font=2, offset=0.1)
 text(0, 0.825, bgc.select, pos=4, cex=0.8, offset=0.1)
 text(-.025, 0.65, "Historical suitabilities", pos=4, cex=0.8, font=2, offset=0.1)
 text(0, 0.525, std.ref, pos=4, cex=0.8, offset=0.1)
 text(-.025, 0.35, paste("Mean projected suitabilities", sep=""), pos=4, cex=0.8, font=2, offset=0.1)
 text(0, 0.225, std.proj, pos=4, cex=0.8, offset=0.1)
 text(-.025, 0.05, paste("Mean rel. suit. richness: ", round(get(metric)[pt]*100), "%", sep=""), pos=4, cex=0.8, offset=0.1, font=2)
 text(0.57, .95, paste("Projected bioclimates"), pos=4, cex=0.8, font=2, offset=0.1)
 
 par(mar=c(3,0,1,0), plt = plt.callout+c(0.225, -0.005, 0.08, -0.04), new = TRUE, mgp=c(0.5,0.1,0))
 bar.bgc <- names(rev(sort(table(BGC.list))))
 bar.zone <- rep(NA, length(bar.bgc))
 for(i in BGCcolors$zone){ bar.zone[grep(i,bar.bgc)] <- i }
 barplot(rev(sort(table(BGC.list))), col=as.character(BGCcolors$HEX[match(bar.zone, BGCcolors$zone)]), horiz=F, las=2, ylab=list("# GCMs", cex=0.7), yaxt="n", cex.names=0.7)
 par(mgp=c(1,0,0))
 axis(2, at=seq(0,12,2), labels=seq(0,12,2), las=1, tck=0, cex.axis=0.8, lwd=0)
 
 ## Callout box for Location 3
 i=3
 bgc.select <- focal.bgcs[i]
 focal.bgc <- which(BGC==focal.bgcs[i])
 pt <- pts.select[i] 

 #convert fractional suitability to standard suitability
 comm.ref.pt <- comm.ref[pt,]
 comm.ref.pt <- round((1-comm.ref.pt)*4+1)
 comm.ref.pt[comm.ref.pt==4] <- 5
 
 #calculate ensemble mean fractional suitability and convert to standard suitability
 for(GCM in GCMs){
   if(GCM==GCMs[1]) comm.proj.pt <- get(paste("comm.proj", GCM, sep="."))[pt,] else comm.proj.pt <- rbind(comm.proj.pt,get(paste("comm.proj", GCM, sep="."))[pt,])
 }
 comm.proj.pt.ensMean <- apply(comm.proj.pt, 2, FUN=mean, na.rm=T)
 comm.proj.pt.ensMean <- round((1-comm.proj.pt.ensMean)*4+1)
 comm.proj.pt.ensMean[comm.proj.pt.ensMean==4] <- 5
 
 std.ref <- paste(paste(names(comm.ref.pt)[which(comm.ref.pt==1)],collapse=""), "(", 
                  paste(names(comm.ref.pt)[which(comm.ref.pt==2)],collapse=""), ")((", 
                  paste(names(comm.ref.pt)[which(comm.ref.pt==3)],collapse=""), "))", sep="")
 std.proj <- paste(paste(names(comm.proj.pt.ensMean)[which(comm.proj.pt.ensMean==1)],collapse=""), "(", 
                   paste(names(comm.proj.pt.ensMean)[which(comm.proj.pt.ensMean==2)],collapse=""), ")((", 
                   paste(names(comm.proj.pt.ensMean)[which(comm.proj.pt.ensMean==3)],collapse=""), "))", sep="")
 
 BGC.list <- vector()
 for(GCM in GCMs){
   BGC.list[which(GCMs==GCM)] <- get(paste("BGC.pred", GCM, rcp, proj.year, sep="."))[pt]
   # print(BGC.pred[pt])
   # print(GCM)
 }
 BGC.list 
 
 plt.callout <- c(0.65, 0.99, 0.81, 0.99)
 
 par(mar=c(0,0,0,0), plt = plt.callout, new = TRUE, mgp=c(2,0.1,0))
 plot(0, xlim=c(0,1), ylim=c(0,1), col="white", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
 rect(-9,-9,9,9, col="white", lwd=0)
 mtext("(C)", side=3, line=-1, adj=-0.08, cex=1, font=2)
 box()
 
 par(mar=c(0,0,0,0), plt = plt.callout, new = TRUE, mgp=c(2,0.1,0))
 plot(0, xlim=c(0,1), ylim=c(0,1), col="white", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
 text(-0.025, .95, paste("Historical bioclimate"), pos=4, cex=0.8, font=2, offset=0.1)
 text(0, 0.825, bgc.select, pos=4, cex=0.8, offset=0.1)
 text(-.025, 0.65, "Historical suitabilities", pos=4, cex=0.8, font=2, offset=0.1)
 text(0, 0.525, std.ref, pos=4, cex=0.8, offset=0.1)
 text(-.025, 0.35, paste("Mean projected suitabilities", sep=""), pos=4, cex=0.8, font=2, offset=0.1)
 text(0, 0.225, std.proj, pos=4, cex=0.8, offset=0.1)
 text(-.025, 0.05, paste("Mean rel. suit. richness: ", round(get(metric)[pt]*100), "%", sep=""), pos=4, cex=0.8, offset=0.1, font=2)
 text(0.57, .95, paste("Projected bioclimates"), pos=4, cex=0.8, font=2, offset=0.1)
 
 par(mar=c(3,0,1,0), plt = plt.callout+c(0.225, -0.005, 0.08, -0.04), new = TRUE, mgp=c(0.5,0.1,0))
 bar.bgc <- names(rev(sort(table(BGC.list))))
 bar.zone <- rep(NA, length(bar.bgc))
 for(i in BGCcolors$zone){ bar.zone[grep(i,bar.bgc)] <- i }
 barplot(rev(sort(table(BGC.list))), col=as.character(BGCcolors$HEX[match(bar.zone, BGCcolors$zone)]), horiz=F, las=2, ylab=list("# GCMs", cex=0.7), yaxt="n", cex.names=0.7)
 par(mgp=c(1,0,0))
 axis(2, at=seq(0,12,2), labels=seq(0,12,2), las=1, tck=0, cex.axis=0.8, lwd=0)
 
 ## Callout box for Location 4
 i=4
 bgc.select <- focal.bgcs[i]
 focal.bgc <- which(BGC==focal.bgcs[i])
 pt <- pts.select[i] 

 #convert fractional suitability to standard suitability
 comm.ref.pt <- comm.ref[pt,]
 comm.ref.pt <- round((1-comm.ref.pt)*4+1)
 comm.ref.pt[comm.ref.pt==4] <- 5
 
 #calculate ensemble mean fractional suitability and convert to standard suitability
 for(GCM in GCMs){
   if(GCM==GCMs[1]) comm.proj.pt <- get(paste("comm.proj", GCM, sep="."))[pt,] else comm.proj.pt <- rbind(comm.proj.pt,get(paste("comm.proj", GCM, sep="."))[pt,])
 }
 comm.proj.pt.ensMean <- apply(comm.proj.pt, 2, FUN=mean, na.rm=T)
 comm.proj.pt.ensMean <- round((1-comm.proj.pt.ensMean)*4+1)
 comm.proj.pt.ensMean[comm.proj.pt.ensMean==4] <- 5
 
 std.ref <- paste(paste(names(comm.ref.pt)[which(comm.ref.pt==1)],collapse=""), "(", 
                  paste(names(comm.ref.pt)[which(comm.ref.pt==2)],collapse=""), ")((", 
                  paste(names(comm.ref.pt)[which(comm.ref.pt==3)],collapse=""), "))", sep="")
 std.proj <- paste(paste(names(comm.proj.pt.ensMean)[which(comm.proj.pt.ensMean==1)],collapse=""), "(", 
                   paste(names(comm.proj.pt.ensMean)[which(comm.proj.pt.ensMean==2)],collapse=""), ")((", 
                   paste(names(comm.proj.pt.ensMean)[which(comm.proj.pt.ensMean==3)],collapse=""), "))", sep="")
 
 BGC.list <- vector()
 for(GCM in GCMs){
   BGC.list[which(GCMs==GCM)] <- get(paste("BGC.pred", GCM, rcp, proj.year, sep="."))[pt]
   # print(BGC.pred[pt])
   # print(GCM)
 }
 BGC.list 
 
 plt.callout <- c(0.01, 0.32, 0.57, 0.75)
 
 par(mar=c(0,0,0,0), plt = plt.callout, new = TRUE, mgp=c(2,0.1,0))
 plot(0, xlim=c(0,1), ylim=c(0,1), col="white", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
 rect(-9,-9,9,9, col="white", lwd=0)
 mtext("(B)", side=3, line=0.25, adj=0.025, cex=1, font=2)
 box()
 
 par(mar=c(0,0,0,0), plt = plt.callout, new = TRUE, mgp=c(2,0.1,0))
 plot(0, xlim=c(0,1), ylim=c(0,1), col="white", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
 text(-0.025, .95, paste("Historical bioclimate"), pos=4, cex=0.8, font=2, offset=0.1)
 text(0, 0.825, bgc.select, pos=4, cex=0.8, offset=0.1)
 text(-.025, 0.65, "Historical suitabilities", pos=4, cex=0.8, font=2, offset=0.1)
 text(0, 0.525, std.ref, pos=4, cex=0.8, offset=0.1)
 text(-.025, 0.35, paste("Mean projected suitabilities", sep=""), pos=4, cex=0.8, font=2, offset=0.1)
 text(0, 0.225, std.proj, pos=4, cex=0.8, offset=0.1)
 text(-.025, 0.05, paste("Mean rel. suit. richness: ", round(get(metric)[pt]*100), "%", sep=""), pos=4, cex=0.8, offset=0.1, font=2)
 text(0.525, .95, paste("Projected bioclimates"), pos=4, cex=0.8, font=2, offset=0.1)
 
 par(mar=c(3,0,1,0), plt = plt.callout+c(0.2, -0.005, 0.08, -0.04), new = TRUE, mgp=c(0.5,0.1,0))
 bar.bgc <- names(rev(sort(table(BGC.list))))
 bar.zone <- rep(NA, length(bar.bgc))
 for(i in BGCcolors$zone){ bar.zone[grep(i,bar.bgc)] <- i }
 barplot(rev(sort(table(BGC.list))), col=as.character(BGCcolors$HEX[match(bar.zone, BGCcolors$zone)]), horiz=F, las=2, ylab=list("# GCMs", cex=0.7), yaxt="n", cex.names=0.7)
 par(mgp=c(1,0,0))
 axis(2, at=seq(0,12,2), labels=seq(0,12,2), las=1, tck=0, cex.axis=0.8, lwd=0)
 
 ## Summary by zone
 par(mar=c(4.5,2,0.1,0.1), plt = c(0.075, 0.365, 0.125, 0.385), new = TRUE, mgp=c(2.5,0.25,0))
 xlim=c(1, length(levels(droplevels(zone[-grep("BAFA|CMA|IMA",zone)]))))
 z <- boxplot(y~zone, ylim=ylim, ylab="", vertical = TRUE, plot=F)
 for(i in 1:length(levels(zone))){ 
   temp <- y[which(zone==levels(zone)[i])]
   z$stats[c(1,5), i] <- quantile(temp[!is.na(temp)],c(0.05, 0.95))
 }
 bxp(z, ylim=ylim, xlim=xlim, xaxt="n", yaxt="n", xaxs="i", ylab="", pch=0,outline=FALSE)
 lines(c(-99,99), c(0,0), lwd=2, col="darkgrey")
 bxp(z, add=T, boxfill = as.character(BGCcolors$HEX[match(levels(droplevels(zone[-grep("BAFA|CMA|IMA",zone)])), BGCcolors$zone)]), xaxt="n", yaxt="n", xaxs="i", ylab="Relative suit. richness", pch=0,outline=FALSE)
 axis(1, at=1:length(levels(zone)), levels(zone), tick=F, las=2, cex=0.8)
 axis(2,at=seq(ylim[1], ylim[2]), labels=paste(round(2^(seq(ylim[1], ylim[2]))*100),"%", sep=""), las=2, tck=0)
 # lines(range(grep("CWH|MH|CDF", levels(zone)))+c(-0.6,0.6), if(scenario=="6variable") rep(2.7,2) else rep(3.2,2), col=alpha("red", 0.5), lwd=3, lty=1)
 # lines(range(grep("ESSF|MS|IDF|PP|BG|ICH|SBPS|SBS|BWBS|SWB", levels(zone)))+c(-0.6,0.6), if(scenario=="6variable") rep(1.5,2) else rep(1.6,2), col=alpha("red", 0.5), lwd=3, lty=1)
 axis(1, at=c(mean(grep("CWH|MH|CDF", levels(zone))), mean(grep("ESSF|MS|IDF|PP|BG|ICH|SBPS|SBS|BWBS|SWB", levels(zone))),mean(grep("CMA|IMA|BAFA", levels(zone)))), 
      labels = c("Coast", "Interior", "Alpine"), tick=F, line=2.95 )
 par(xpd=NA)
 x1 <- range(grep("CWH|MH|CDF", levels(zone)))+c(-0.4,0.4)
 x2 <- range(grep("ESSF|MS|IDF|PP|BG|ICH|SBPS|SBS|BWBS|SWB", levels(zone)))+c(-0.4,0.4)
 # x3 <- range(grep("CMA|IMA|BAFA", levels(zone)))+c(-0.4,0.4)
 bracketpos=-5.3
 lines(rep(x1, each=2), bracketpos*c(0.95, 1,1,0.95))
 lines(rep(x2, each=2), bracketpos*c(0.95, 1,1,0.95))
 # lines(rep(x3, each=2), bracketpos*c(0.95, 1,1,0.95))
 par(xpd=F)  
 mtext("(F)", side=3, line=-1.5, adj=0.025, cex=1, font=2)
 
 
 dev.off()
 
 print(edatope)
}
  print(proj.year)
}

 
 # # #####################
 # # #### exploratory Suitability profile for individual locations
 # # #####################
 # 
 # edatope=edatopes[2]
 # rcp=rcps[1]
 # proj.year=proj.years[2]
 # 
 # 
 # #select a single gridpoint within a specified BGC
 # bgc.select <- "MSdk"
 # focal.bgc <- which(BGC==bgc.select)
 # par(mar=c(4,4,1,1))
 # hist(SuitTurnover[focal.bgc])
 # # hist(SuitPersistence[focal.bgc])
 # i=0.985
 # for(i in seq(0.805, 0.995, 0.01)){
 #   q.select <- i #quantile of turnover to select the example grid point from 
 #   pt <- focal.bgc[which(SuitTurnover[focal.bgc]==min(SuitTurnover[focal.bgc][SuitTurnover[focal.bgc]>quantile(SuitTurnover[focal.bgc], q.select, na.rm=T)], na.rm=T))[1]] #select the grid point
 #   pt <- pt[1]
 #   # pt <- plotOrder[!is.na(plotOrder)][161020]
 #   #convert fractional suitability to standard suitability
 #   comm.ref.pt <- comm.ref[pt,]
 #   comm.ref.pt <- round((1-comm.ref.pt)*4+1)
 #   comm.ref.pt[comm.ref.pt==4] <- 5
 #   
 #   #calculate ensemble mean fractional suitability and convert to standard suitability
 #   for(GCM in GCMs){
 #     if(GCM==GCMs[1]) comm.proj.pt <- get(paste("comm.proj", GCM, sep="."))[pt,] else comm.proj.pt <- rbind(comm.proj.pt,get(paste("comm.proj", GCM, sep="."))[pt,])
 #   }
 #   comm.proj.pt.ensMean <- apply(comm.proj.pt, 2, FUN=mean, na.rm=T)
 #   comm.proj.pt.ensMean <- round((1-comm.proj.pt.ensMean)*4+1)
 #   comm.proj.pt.ensMean[comm.proj.pt.ensMean==4] <- 5
 #   
 #   std.ref <- paste(paste(names(comm.ref.pt)[which(comm.ref.pt==1)],collapse=""), "(", 
 #                    paste(names(comm.ref.pt)[which(comm.ref.pt==2)],collapse=""), ")((", 
 #                    paste(names(comm.ref.pt)[which(comm.ref.pt==3)],collapse=""), "))", sep="")
 #   std.proj <- paste(paste(names(comm.proj.pt.ensMean)[which(comm.proj.pt.ensMean==1)],collapse=""), "(", 
 #                     paste(names(comm.proj.pt.ensMean)[which(comm.proj.pt.ensMean==2)],collapse=""), ")((", 
 #                     paste(names(comm.proj.pt.ensMean)[which(comm.proj.pt.ensMean==3)],collapse=""), "))", sep="")
 #   
 #   GCM <- GCMs[2]
 #   pt.bgc <- get(paste("BGC.pred", GCM, rcp, proj.year, sep="."))[pt]
 #   
 #   #map
 #   par(mar=c(0,0,0,0))
 #   values(X) <- SuitTurnover[plotOrder]
 #   
 #   # breakpoints <- c(-99,seq(0,2,0.25),99); length(breakpoints)
 #   # ColScheme <- brewer.pal(10,"RdBu"); length(ColScheme)
 #   breakpoints <- c(seq(0,1,.1),100); length(breakpoints)
 #   ColScheme <- c(brewer.pal(9,"YlGnBu"), "black", "red"); length(ColScheme)
 #   
 #   std.ref
 #   std.proj
 #   # quantile(SuitTurnover[focal.bgc], q.select)
 #   # SuitPersistence[pt]
 #   # SuitTurnover[pt]
 #   
 #   BGC.list <- vector()
 #   for(GCM in GCMs){
 #     BGC.list[which(GCMs==GCM)] <- get(paste("BGC.pred", GCM, rcp, proj.year, sep="."))[pt]
 #     # print(BGC.pred[pt])
 #     # print(GCM)
 #   }
 #   BGC.list 
 #   
 #   
 #   par(mfrow=c(1,1), mar=c(0,0,0,0))
 #   image(X, col=ColScheme, breaks=breakpoints, maxpixels= ncell(X), ylim=extent(X)[3:4]-c(200000, 0))
 #   # plot(bdy.bc, add=T, border="black", lwd=0.4)
 #   # points(pts.df[plotOrder[!is.na(plotOrder)][161025],], cex=2, lwd=2)
 #   points(pts.df[pt,], cex=2, lwd=2)
 #   legend("topright", legend=i, cex=2)
 #   
 #   par(mar=c(0,0,0,0), plt = c(0.01, 0.335, 0.01, 0.38), new = TRUE, mgp=c(2,0.1,0))
 #   plot(0, xlim=c(0,1), ylim=c(0,1), col="white", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
 #   box()
 #   
 #   par(mar=c(0,0,0,0), plt = c(0.01, 0.335, 0.01, 0.38), new = TRUE, mgp=c(2,0.1,0))
 #   plot(0, xlim=c(0,1), ylim=c(0,1), col="white", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
 #   text(-0.025, .95, paste("Historical bioclimate"), pos=4, cex=0.8, font=2, offset=0.1)
 #   text(0, 0.85, bgc.select, pos=4, cex=0.8, offset=0.1)
 #   text(-.025, 0.75, "Historical suitabilities", pos=4, cex=0.8, font=2, offset=0.1)
 #   text(0, 0.65, std.ref, pos=4, cex=0.8, offset=0.1)
 #   text(-.025, 0.55, paste("Mean projected suitabilities", sep=""), pos=4, cex=0.8, font=2, offset=0.1)
 #   text(0, 0.45, std.proj, pos=4, cex=0.8, offset=0.1)
 #   text(-.025, 0.35, paste("Mean turnover: ", round(SuitTurnover[pt]*100), "%", sep=""), pos=4, cex=0.8, offset=0.1, font=2)
 #   text(-.025, 0.25, paste("Mean persistence: ", round(SuitPersistence[pt]*100), "%", sep=""), pos=4, cex=0.8, offset=0.1, font=2)
 #   text(-.025, 0.15, paste("Mean richness change: ", round(SuitRichnessChangePct[pt]*100), "%", sep=""), pos=4, cex=0.8, offset=0.1, font=2)
 #   
 #   par(mar=c(3,0,1,0), plt = c(0.8, 0.99, 0.09, 0.15), new = TRUE, mgp=c(0.5,0.1,0))
 #   barplot(rev(sort(table(BGC.list))), horiz=F, las=2, ylab=list("# of GCMs", cex=0.7), yaxt="n", cex.names=0.6)
 #   par(mgp=c(1,0,0))
 #   axis(2, las=1, tck=0, cex.axis=0.7)
 #   title(main=list("Projected bioclimates", cex=0.7))
 # }
 # # SuitPersistence examples
 # # decline: IDFxh4, (CWHmm2), CWHws1
 # # stable: 
 # # increase: ESSFdkp
 # 
 # 
 # # SuitTurnover examples
 # # Low: SBSmk1
 # # High: ESSFwc4 (includes some CWHws1/2 and MHmm2 units); 
 # # High: ICHxw; 
 # 
 # levels(bymedian)
 # 
 # 