
##======================================================================================
## CCISS Publication Scripts
## Step 4d - Figures - Suitability for individual species and groups of species
##======================================================================================


source("./_CCISS_Packages.R") ## packages required
source("./_CCISS_Functions.R") ## common functions
source("./_CCISS_Parameters.R") ## settings used through all scripts

#===============================================================================
# create a dem from the climateBC input data
#===============================================================================

## create a dem from the climateBC input data
points <- read.csv(paste("inputs\\",grid,".csv", sep=""))
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
bdy.bc <- readOGR("inputs\\shapes\\ProvincialOutline.shp")


#===============================================================================
# Import suitability tables
#===============================================================================

S1 <- treesuit
S1 <- unique(S1[,1:4])

# select the species to run the analysis on
spps <- unique(S1$Spp)
spps <- spps[-which(spps=="X")]
spps.candidate <- spps.lookup$TreeCode[-which(spps.lookup$Exclude=="x")]
spps <- spps[which(spps%in%spps.candidate)] 

## non-THLB BGCs for exclusion from results
BGC <- points$ID2
BGC <- gsub(" ","",BGC)
BGCs_notin_THLB <- BGCs_notin_THLB$BGC[which(BGCs_notin_THLB$Exlude=="x")]
exclude <- which(BGC%in%BGCs_notin_THLB)

#BGC zones
zone <- rep(NA, length(BGC))
for(i in BGCcolors.BC$zone){ zone[grep(i,BGC)] <- i }
table(zone)
zone <- factor(zone, levels=c("CDF", "CWH", "MH", "ESSF", "MS", "IDF", "PP", "BG", "ICH", "SBPS", "SBS", "BWBS", "SWB", "CMA", "IMA", "BAFA"))


#===============================================================================
# Manuscript plot for a subset of species
#===============================================================================
rcp=rcps[1]
hist.year=hist.years[2]
proj.year=proj.years[2]

edatope="C4"
# for(proj.year in proj.years[1:2]){
# for(edatope in edatopes){

spps.matrix <- matrix(c("Pl", "Fd", "Cw","Ba", "Bl", "Bg", "Yc", "Ss", "Hm", "Lw", "Hw", "Py", "Dr", "Ep", "At"), 5, byrow=T)  

# x11(width=6.5, height=8.5, pointsize=12)

matrow=1
for(matrow in 1:dim(spps.matrix)){
 spps <- spps.matrix[matrow,]

# png(filename=paste("results\\Manu_Suitability_Groups\\CCISS.manu.Suitability",spps[1],spps[2],spps[3], edatope, rcp, proj.year,"png",sep="."), type="cairo", units="in", width=6.5, height=8.5, pointsize=12, res=600)
pdf(file=paste("results\\CCISS.Fig3.Suitability",spps[1],spps[2],spps[3], edatope, rcp, proj.year,"pdf",sep="."), width=6.5, height=8.5, pointsize=12)
par(mar=c(0,0,0,0), mfrow=c(3,1), bg="white")

for(spp in spps){

  RefSuit <- read.csv(paste("outputs\\Suit.ref", grid, spp, edatope, "csv", sep="."))[,1]
  outRange.base <- RefSuit==5
  RefSuit[RefSuit==5] <- 4
  RefSuit[is.na(RefSuit)] <- 4
  table(RefSuit)
  
  HistSuit <- read.csv(paste("outputs\\Suit", grid, hist.year, spp, edatope, "csv", sep="."))[,1]
  HistSuit[HistSuit==5] <- 4
  HistSuit[is.na(HistSuit)] <- 4
  ChangeSuit.hist <- RefSuit-HistSuit
  
    # compile the GCM projections into a data frame
      ProjSuit <- data.frame(temp=rep(NA, length(RefSuit))) #initiate the data frame with a dummy column
      ChangeSuit <- data.frame(temp=rep(NA, length(RefSuit))) #initiate the data frame with a dummy column
      for(GCM in GCMs){
        temp <- read.csv(paste("outputs\\Suit", grid, GCM, rcp, proj.year, spp, edatope, "csv", sep="."))
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
  panel <- paste("(", letters[which(spps==spp)],")", sep="")
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
  text(xl-30000, mean(c(yb,yt))-30000, paste("Change to feasible/unfeasible\n(", proj.year.name[which(proj.years==proj.year)], "), % of GCMs", sep=""), srt=90, pos=3, cex=0.85, font=2)
  mtext(paste("(", letters[c(2,6,10)][which(spps==spp)],")", sep=""), side=3, line=-2.5, adj=0.22, cex=0.8, font=2)
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
         fill=ColScheme, bty="n", cex=0.9, title="Historical feasibility", inset=0.015)
  # }
    mtext(paste("(", letters[c(1,5,9)][which(spps==spp)],")", sep=""), side=3, line=-3.75, adj=0.05, cex=0.8, font=2)
    
    # ##=================================
    # # recent period binary change
    # binary <- rep(0, length(RefSuit))
    # binary[outRange.base==T] <- NA
    # binary[outRange.base==T][HistSuit[outRange.base==T] < 4] <- 1   
    # binary[outRange.base==F][HistSuit[outRange.base==F] == 4] <- -1
    # 
    # values(X) <- binary[plotOrder]
    # 
    # ColScheme <- c(brewer.pal(11,"RdBu")[2], "grey90", brewer.pal(11,"RdBu")[10]); length(ColScheme)
    # 
    # par(plt = c(0.6, 0.95, 0.3, 1), new = TRUE)
    # plot(bdy.bc, border="black", lwd=0.4)
    # image(X, add=T, xaxt="n", yaxt="n", col=ColScheme, maxpixels= ncell(X))
    # legend(1400000, 1550000, legend=c("Expand", "Persist", "Retreat"), 
    #        fill=rev(ColScheme), bty="n", cex=0.9, title=paste("Observed\nclimate\n(",hist.year.name[which(hist.years==hist.year)],")", sep=""), inset=0.015)
    # mtext(paste("(", letters[c(3,7,11)][which(spps==spp)],")", sep=""), side=3, line=-3.25, adj=0.1, cex=0.8, font=2)
    
  
    ##=================================
    # ALTERNATE: map of suitability change
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
    text(xl-30000, mean(c(yb,yt))-30000, paste("Mean change\nin feasibility (", proj.year.name[which(proj.years==proj.year)], ")", sep=""), srt=90, pos=3, cex=0.9, font=2)
    # }
    par(xpd=F)
    mtext(paste("(", letters[c(3,7,11)][which(spps==spp)],")", sep=""), side=3, line=-3.25, adj=0.1, cex=0.8, font=2)

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
    bxp(z, add=T, boxfill = as.character(BGCcolors.BC$HEX[match(levels(zone), BGCcolors.BC$zone)]), xaxt="n", yaxt="n", xaxs="i", ylab="", pch=0,outline=FALSE)
    axis(1, at=1:length(levels(zone)), levels(zone), tick=F, las=2, cex.axis=0.8)
    axis(2,at=seq(ylim[1], ylim[2], 3), seq(ylim[1], ylim[2], 3), las=2, tck=0)
    mtext(paste("(", letters[c(4,8,12)][which(spps==spp)],")", sep=""), side=3, line=1, adj=0.975, cex=0.8, font=2)
    mtext(paste("Mean feasibility change (", proj.year.name[which(proj.years==proj.year)], ")", sep=""), side=3, line=0.1, adj=.975, cex=0.5, font=2)

    print(spp)
}
dev.off()
print(matrow)
}

#===============================================================================
# Manuscript plot of all edatopes for one species
#===============================================================================

# select the species to run the analysis on
spps <- unique(S1$Spp)
spps <- spps[-which(spps=="X")]
spps.candidate <- spps.lookup$TreeCode[-which(spps.lookup$Exclude=="x")]
spps <- spps[which(spps%in%spps.candidate)] 

spp="Sx"

for(proj.year in proj.years){
  for(spp in spps){
    
  # png(filename=paste("results\\Manu_Suitability_Indiv\\CCISS.manu.Suitability.edatopes",spp, rcp, proj.year,"png",sep="."), type="cairo", units="in", width=6.5, height=8.5, pointsize=12, res=600)
  pdf(file=paste("results\\Manu_Suitability_Indiv\\CCISS.Suitability_Indiv",spp, rcp, proj.year,"pdf",sep="."), width=6.5, height=8.5, pointsize=12)
  par(mar=c(0,0,0,0), mfrow=c(3,1), bg="white")
  
  for(edatope in edatopes){
    
    RefSuit <- read.csv(paste("outputs\\Suit.ref", grid, spp, edatope, "csv", sep="."))[,1]
    outRange.base <- RefSuit==5
    RefSuit[RefSuit==5] <- 4
    RefSuit[is.na(RefSuit)] <- 4
    
    HistSuit <- read.csv(paste("outputs\\Suit", grid, hist.year, spp, edatope, "csv", sep="."))[,1]
    HistSuit[HistSuit==5] <- 4
    HistSuit[is.na(HistSuit)] <- 4
    ChangeSuit.hist <- RefSuit-HistSuit
    
    # compile the GCM projections into a data frame
    ProjSuit <- data.frame(temp=rep(NA, length(RefSuit))) #initiate the data frame with a dummy column
    ChangeSuit <- data.frame(temp=rep(NA, length(RefSuit))) #initiate the data frame with a dummy column
    for(GCM in GCMs){
      temp <- read.csv(paste("outputs\\Suit", grid, GCM, rcp, proj.year, spp, edatope, "csv", sep="."))
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
    panel <- paste("(", letters[which(spps==spp)],")", sep="")
    mtext(paste("Site type: ", edatope, " (", edatope.name[which(edatopes==edatope)], ")", sep=""), side=3, line=-1.55, adj=0.01, cex=0.75, font=2)
    mtext(if(spp%in%spps.lookup$TreeCode) bquote(bold(.(spp))~"-"~.(Common)) else bquote(bold(.(spp))), side=3, line=-2.75, adj=0.01, cex=0.7, font=1)
    box()
    
    ##=================================
    # FUTURE MAP of binary appearance/disappearance
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
    text(xl-30000, mean(c(yb,yt))-30000, paste("Change to feasible/unfeasible\n(", proj.year.name[which(proj.years==proj.year)], ") % of GCMs", sep=""), srt=90, pos=3, cex=0.9, font=2)
    mtext(paste("(", letters[c(2,6,10)][which(edatopes==edatope)],")", sep=""), side=3, line=-2.75, adj=0.22, cex=0.8, font=2)
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
           fill=ColScheme, bty="n", cex=0.9, title="Historical feasibility", inset=0.015)
    # }
    mtext(paste("(", letters[c(1,5,9)][which(edatopes==edatope)],")", sep=""), side=3, line=-3.75, adj=0.05, cex=0.8, font=2)
    
    # ##=================================
    # # RECENT PERIOD binary change
    # binary <- rep(0, length(RefSuit))
    # binary[outRange.base==T] <- NA
    # binary[outRange.base==T][HistSuit[outRange.base==T] < 4] <- 1   
    # binary[outRange.base==F][HistSuit[outRange.base==F] == 4] <- -1
    # 
    # values(X) <- binary[plotOrder]
    # 
    # ColScheme <- c(brewer.pal(11,"RdBu")[2], "grey90", brewer.pal(11,"RdBu")[10]); length(ColScheme)
    # 
    # par(plt = c(0.6, 0.95, 0.3, 1), new = TRUE)
    # plot(bdy.bc, border="black", lwd=0.4)
    # image(X, add=T, xaxt="n", yaxt="n", col=ColScheme, maxpixels= ncell(X))
    # legend(1400000, 1550000, legend=c("Expand", "Persist", "Retreat"), 
    #        fill=rev(ColScheme), bty="n", cex=0.9, title=paste("Observed\nclimate\n(",hist.year.name[which(hist.years==hist.year)],")", sep=""), inset=0.015)
    # mtext(paste("(", letters[c(3,7,11)][which(edatopes==edatope)],")", sep=""), side=3, line=-3.25, adj=0.1, cex=0.8, font=2)

    ##=================================
    # ALTERNATE: map of suitability change
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
    text(xl-30000, mean(c(yb,yt))-30000, paste("Mean change\nin feasibility (", proj.year.name[which(proj.years==proj.year)], ")", sep=""), srt=90, pos=3, cex=0.9, font=2)
    # }
    par(xpd=F)
    mtext(paste("(", letters[c(3,7,11)][which(edatopes==edatope)],")", sep=""), side=3, line=-3.25, adj=0.1, cex=0.8, font=2)
    
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
    bxp(z, add=T, boxfill = as.character(BGCcolors.BC$HEX[match(levels(zone), BGCcolors.BC$zone)]), xaxt="n", yaxt="n", xaxs="i", ylab="", pch=0,outline=FALSE)
    axis(1, at=1:length(levels(zone)), levels(zone), tick=F, las=2, cex.axis=0.8)
    axis(2,at=seq(ylim[1], ylim[2], 3), seq(ylim[1], ylim[2], 3), las=2, tck=0)
    mtext(paste("(", letters[c(4,8,12)][which(edatopes==edatope)],")", sep=""), side=3, line=1, adj=0.975, cex=0.8, font=2)
    mtext(paste("Mean feasibility change (", proj.year.name[which(proj.years==proj.year)], ")", sep=""), side=3, line=0.1, adj=.975, cex=0.55, font=2)
    
    print(edatope)
  }
  dev.off()
  print(spp)
}
print(proj.year)
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
for(proj.year in proj.years[-3]){
  for(edatope in edatopes){
    
    for(spp in spps[c(2,3,6,7,12,14,24)]){
      
      png(filename=paste("results\\SI_Suitability_Basic\\CCISS.manu.Suitability.basic",spp, edatope, rcp, proj.year,"png",sep="."), type="cairo", units="in", width=6.5, height=6, pointsize=10, res=600)
      par(mar=c(0,0,0,0), mfrow=c(2,2), bg="white")
      
      RefSuit <- read.csv(paste("outputs\\Suit.ref", grid, spp, edatope, "csv", sep="."))[,1]
      outRange.base <- RefSuit==5
      RefSuit[RefSuit==5] <- 4
      RefSuit[is.na(RefSuit)] <- 4
      
      # compile the GCM projections into a data frame
      ProjSuit <- data.frame(temp=rep(NA, length(RefSuit))) #initiate the data frame with a dummy column
      ChangeSuit <- data.frame(temp=rep(NA, length(RefSuit))) #initiate the data frame with a dummy column
      for(GCM in GCMs){
        temp <- read.csv(paste("outputs\\Suit", grid, GCM, rcp, proj.year, spp, edatope, "csv", sep="."))
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
             fill=ColScheme, bty="n", cex=0.9, title="Historical feasibility", inset=0.015)
      # }
      # box()
      Common <- as.character(spps.lookup$EnglishName[which(spps.lookup$TreeCode==spp)])
      Latin <- as.character(spps.lookup$ScientificName[which(spps.lookup$TreeCode==spp)])
      panel <- paste("(", letters[which(spps==spp)],")", sep="")
      mtext(if(spp%in%spps.lookup$TreeCode) bquote(bold(.(spp))~"-"~.(Common)) else bquote(.(panel)~bold(.(spp))),
            side=1, line=-1.75, adj=0.01, cex=0.8, font=2)
      # mtext(if(spp%in%spps.lookup$TreeCode) bquote(.(panel)~bold(.(spp))~"-"~.(Common)~"("*italic(.(Latin)*")")) else bquote(.(panel)~bold(.(spp))),
      #       side=3, line=-1.75, adj=0.01, cex=0.8, font=2)
      mtext(paste("Site type: ", edatope, " (", edatope.name[which(edatopes==edatope)], ")", sep=""), side=1, line=-2.75, adj=0.01, cex=0.7, font=1)
      
      # map of suitability change
      breakpoints <- seq(-3,3,0.5); length(breakpoints)
      labels <- c("-3","-2", "-1", "no\nchange", "+1","+2","+3")
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
      text(xl-30000, mean(c(yb,yt))-30000, paste("Mean change\nin feasibility (", proj.year.name[which(proj.years==proj.year)], ")", sep=""), srt=90, pos=3, cex=0.9, font=2)
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
      labels <- c("Low", "Medium", "High")
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
      text(xl-30000, mean(c(yb,yt))-30000, paste("Ensemble disagreement\non trend"), srt=90, pos=3, cex=1, font=2)
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



#===============================================================================
# presentation plot
#===============================================================================
rcp=rcps[1]
hist.year=hist.years[2]
proj.year=proj.years[2]

edatope="C4"

# x11(width=6.5, height=8.5, pointsize=12)

# for(proj.year in proj.years[-3]){
  for(edatope in edatopes){
    # spp="Sx"
    for(spp in spps){
      
  png(filename=paste("results\\SI_Suitability_Presentation\\CCISS.manu.Suitability",spp, edatope, rcp, proj.year,"png",sep="."), type="cairo", units="in", width=7.75, height=5.5, pointsize=12, res=600)
  # pdf(file=paste("results\\CCISS.Fig3.Suitability",spps[1],spps[2],spps[3], edatope, rcp, proj.year,"pdf",sep="."), width=6.5, height=8.5, pointsize=12)
  par(mar=c(0,0,0,0), mfrow=c(1,1), bg="white")
  
    RefSuit <- read.csv(paste("outputs\\Suit.ref", grid, spp, edatope, "csv", sep="."))[,1]
    outRange.base <- RefSuit==5
    RefSuit[RefSuit==5] <- 4
    RefSuit[is.na(RefSuit)] <- 4
    table(RefSuit)
    
    HistSuit <- read.csv(paste("outputs\\Suit", grid, hist.year, spp, edatope, "csv", sep="."))[,1]
    HistSuit[HistSuit==5] <- 4
    HistSuit[is.na(HistSuit)] <- 4
    ChangeSuit.hist <- RefSuit-HistSuit
    
    # compile the GCM projections into a data frame
    ProjSuit <- data.frame(temp=rep(NA, length(RefSuit))) #initiate the data frame with a dummy column
    ChangeSuit <- data.frame(temp=rep(NA, length(RefSuit))) #initiate the data frame with a dummy column
    for(GCM in GCMs){
      temp <- read.csv(paste("outputs\\Suit", grid, GCM, rcp, proj.year, spp, edatope, "csv", sep="."))
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
    plot(0, col="white", xaxt="n", yaxt="n", xlab="", ylab="", axes=F)
    Common <- as.character(spps.lookup$EnglishName[which(spps.lookup$TreeCode==spp)])
    Latin <- as.character(spps.lookup$ScientificName[which(spps.lookup$TreeCode==spp)])
    panel <- paste("(", letters[which(spps==spp)],")", sep="")
    mtext(if(spp%in%spps.lookup$TreeCode) bquote(bold(.(spp))~"-"~.(Common)) else bquote(bold(.(spp))),
          side=3, line=-5.75, adj=0.01, cex=1.2, font=2)
    # mtext(if(spp%in%spps.lookup$TreeCode) bquote(.(panel)~bold(.(spp))~"-"~.(Common)~"("*italic(.(Latin)*")")) else bquote(.(panel)~bold(.(spp))),
    #       side=3, line=-1.75, adj=0.01, cex=0.8, font=2)
    mtext(paste("Site type: ", edatope, " (", edatope.name[which(edatopes==edatope)], ")", sep=""), side=3, line=-6.75, adj=0.01, cex=1, font=1)
    # box()
    
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
    
    par(plt = c(0.175,1,0,1), new = TRUE)
    plot(bdy.bc, border="black", lwd=0.4)
    image(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks=breakpoints, maxpixels= ncell(X))
    # mtext(paste("Edatope:", edatope), side=1, line=-1.5, adj=0.02, cex=1.1, font=2)
    # if(spp==spps[1]){
    xl <- 1600000; yb <- 1000000; xr <- 1700000; yt <- 1700000
    rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme)
    text(rep(xr+10000,length(labels)),seq(yb,yt,(yt-yb)/(length(GCMs)-1))[c(3,9)],labels,pos=4,cex=0.9,font=0.8, srt=90)
    text(rep(xr-20000,length(labels)),seq(yb,yt,(yt-yb)/(length(GCMs)-1))[c(1,8,15)],c("100%", "0%", "100%"),pos=4,cex=0.8,font=1)
    text(xl-30000, mean(c(yb,yt))-30000, paste("Change to feasible/unfeasible\n(", proj.year.name[which(proj.years==proj.year)], "), % of GCMs", sep=""), srt=90, pos=3, cex=0.99, font=2)
    # mtext(paste("(", letters[c(2,6,10)][which(spps==spp)],")", sep=""), side=3, line=-2.5, adj=0.22, cex=0.8, font=2)
    # legend("bottomleft", legend=c(spp, paste("Edatope:", edatope), proj.year, rcp, " "), cex=1.4, bty="n", inset=-0.05)
    # }
    # box()
    
    ##=================================
    # historical suitability
    breakseq <- c(0.5,1.5,2.5,3.5,5)
    ColScheme <- c(brewer.pal(9,"Greys")[9], brewer.pal(9,"Greens")[7], brewer.pal(9,"Greens")[4], "white")
    ColScheme <- c("darkgreen", "dodgerblue1", "gold2", "white")
    length(ColScheme)
    
    par(plt = c(0, 0.475, 0, 0.55), new = TRUE)
    values(X) <- RefSuit[plotOrder]
    plot(bdy.bc, border="black", lwd=0.4)
    image(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks=breakseq, maxpixels= ncell(X))
    plot(bdy.bc, add=T, border="black", lwd=0.4)
    # if(spp==spps[1]){
    par(xpd=T)
    legend("topleft", legend=c("1 (primary)", "2 (secondary)", "3 (tertiary)"), 
           fill=ColScheme, bty="n", cex=0.9, title="Historical feasibility", inset=c(0,-0.25),)
    par(xpd=F)
    # }
    # mtext(paste("(", letters[c(1,5,9)][which(spps==spp)],")", sep=""), side=3, line=-3.75, adj=0.05, cex=0.8, font=2)
    

  dev.off()
  print(spp)
    }
    print(proj.year)
  # }
  print(edatope)
}
    
