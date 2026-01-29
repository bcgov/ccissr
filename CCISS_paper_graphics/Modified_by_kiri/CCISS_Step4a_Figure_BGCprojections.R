
##======================================================================================
## CCISS Publication Scripts
## Step 4a - Figure - BGC projections for British Columbia
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


setwd("C:\\Colin\\SpatialData\\Boundaries")

###country boundaries
# ORIGIONAL SOURCE: http://www.diva-gis.org/gdata
countries <- readOGR(dsn="countries", layer='countries')
countries.NA <- countries[grep("Canada|United States|Mexico", countries$COUNTRY),]

# ####### create a polygon mask for North America. 
# my_box = as(extent(-179, -50, -20, 84), "SpatialPolygons")      		# convert extent box to shapefile (rectangle)
# cont.NA <- unionSpatialPolygons(countries.NA, rep(1,length(countries.NA$OBJECTID)))
# cont.NA.g <- gSimplify(cont.NA, tol=0.01, topologyPreserve=TRUE)
# proj4string(my_box) = projection(cont.NA)				# assign spatial projection to extent object
# mask.NA <- gDifference(my_box, cont.NA.g)
# projection(mask.NA)  #verify latlong projection of the study area boundary
P4S.latlon <- CRS("+proj=longlat +datum=WGS84")

### admin boundaries
bdy.usa1 <- readOGR("USA_adm",'USA_adm1')
bdy.usa <- gSimplify(bdy.usa1, tol=0.01, topologyPreserve=TRUE) #generalize the linework

# areas <- lapply(bdy.usa@polygons, function(x) sapply(x@Polygons, function(y) y@area))
# bigpolys <- lapply(areas, function(x) which(x > 1000))
# for(i in 1:length(bigpolys)){
#   if(length(bigpolys[[i]]) >= 1 && bigpolys[[i]][1] >= 1){
#     bdy.usa@polygons[[i]]@Polygons <- bdy.usa@polygons[[i]]@Polygons[bigpolys[[i]]]
#     bdy.usa@polygons[[i]]@plotOrder <- 1:length(bdy.usa@polygons[[i]]@Polygons)
#   }
# }

#===============================================================================
# Set analysis Parameters
#===============================================================================

setwd("C:\\Colin\\Projects\\2019_CCISS")

grid <- "BC2kmGrid"

GCMs <-  c("ACCESS1-0","CanESM2","CCSM4","CESM1-CAM5","CNRM-CM5","CSIRO-Mk3-6-0", "GFDL-CM3","GISS-E2R", "HadGEM2-ES", "INM-CM4", "IPSL-CM5A-MR", "MIROC-ESM", "MIROC5", "MPI-ESM-LR","MRI-CGCM3")
rcps <- c("rcp45", "rcp85")
proj.years <- c(2025, 2055, 2085)
hist.years <- c(1995, 2004, 2005, 2009, 2014, 2017)
hist.year.name <- c("1991-2000", "1991-2017", "2001-2010", "2001-2017","2011-2017", "2017")
edatopes<- c("B2", "C4", "D6")
edatope.name <- c("Subxeric-poor", "Mesic-medium", "Hygric-rich")
proj.year.name=c("2020s", "2050s", "2080s")
rcp.name=c("RCP4.5", "RCP8.5")
BGCcolors <- read.csv("C:\\Colin\\Projects\\2019_CCISS\\InputData\\BGCzone_Colorscheme.csv")
rcp.focal="rcp45"
proj.year.focal=2025


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
points <- read.csv(paste("InputData\\",grid,".csv", sep=""))
BGC <- points$ID2
BGC <- gsub(" ","",BGC)  
zone <- rep(NA, length(BGC))
for(i in BGCcolors$zone){ zone[grep(i,BGC)] <- i }

## reference period BGC
BGC.pred.ref <- as.character(read.csv(paste("OutputData\\BGC.pred", grid, "ref.csv", sep="."))[,1])
zone.pred.ref <- rep(NA, length(BGC))
for(i in BGCcolors$zone){ zone.pred.ref[grep(i,BGC.pred.ref)] <- i }

# Historical BGC
for(hist.year in hist.years){
  BGC.pred <- as.character(read.csv(paste("OutputData\\BGC.pred", grid,hist.year,"csv", sep="."))[,1])
  assign(paste("BGC.pred", hist.year, sep="."), BGC.pred) #bgc projection
  
  #did the BGC change?
  assign(paste("BGC.change", hist.year, sep="."), BGC.pred!=BGC) #did the BGC unit change (mapped baseline)
  assign(paste("BGC.change.pred", hist.year, sep="."), BGC.pred!=BGC.pred.ref) #did the BGC unit change (predicted baseline)
  
  #did the zone change? 
  zone.pred <- rep(NA, length(BGC.pred))
  for(i in BGCcolors$zone){ zone.pred[grep(i,BGC.pred)] <- i }
  zone.pred[-which(BGC.pred%in%unique(BGC))] <- "Exotic"
  assign(paste("zone.change", hist.year, sep="."), zone.pred!=zone) #did the BGC zone change (mapped baseline)
  assign(paste("zone.change.pred", hist.year, sep="."), zone.pred!=zone.pred.ref) #did the BGC zone change (predicted baseline)
  
  print(hist.year)
}

# Future BGC
PredSum <- data.frame()
for(rcp in rcps){
  for(proj.year in proj.years){
    for(GCM in GCMs){
      BGC.pred <- as.character(read.csv(paste("OutputData\\BGC.pred",grid, GCM, rcp, proj.year,".csv", sep=""))[,1])
      assign(paste("BGC.pred", GCM, rcp, proj.year, sep="."), BGC.pred) #bgc projection
      PredSum <- rbind(PredSum, as.data.frame(table(BGC.pred)))
      
      #did the BGC change?
      assign(paste("BGC.change", GCM, rcp, proj.year, sep="."), BGC.pred!=BGC) #did the BGC unit change (mapped baseline)
      assign(paste("BGC.change.pred", GCM, rcp, proj.year, sep="."), BGC.pred!=BGC.pred.ref) #did the BGC unit change (predicted baseline)
      
      #did the zone change? 
      zone.pred <- rep(NA, length(BGC.pred))
      for(i in BGCcolors$zone){ zone.pred[grep(i,BGC.pred)] <- i }
      zone.pred[-which(BGC.pred%in%unique(BGC))] <- "Exotic"
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
# for(rcp in rcps){
  proj.year=proj.year.focal
  # for(proj.year in proj.years){
    temp <- as.data.frame(matrix(rep(NA, length(BGC.pred)*length(GCMs)), nrow=length(BGC.pred), ncol=length(GCMs)))
    for(GCM in GCMs){
      BGC.pred <- get(paste("BGC.pred", GCM, rcp, proj.year, sep="."))
      #add votes to votes matrix
      temp[,which(GCMs==GCM)] <- BGC.pred
      print(GCM)
    }
    vote.winner <- function(x){return(names(which(table(x)==max(table(x))))[1])}
    agreement <- function(x){return(max(table(x)))}
    assign(paste("BGC.pred.ensemble", rcp, proj.year, sep="."), apply(temp, 1, vote.winner))
    # assign(paste("BGC.pred.agreement", rcp, proj.year, sep="."), apply(temp, 1, agreement))
    
#     print(proj.year)
#   }
#   print(rcp)
# }

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


#===============================================================================
# create a dem from the climateBC input data
#===============================================================================

    setwd("C:\\Colin\\Projects\\2019_CCISS")
    
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

X <- dem
values(X) <- points$el[df$z]

#===============================================================================
# generate the vector of mapped BGCs
#===============================================================================

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
# generic spatial data
#===============================================================================
### admin boundaries
bdy.bc <- readOGR("InputData\\BC_AB_US_Shp\\ProvincialOutline.shp")
bdy.us <- readOGR("InputData\\BC_AB_US_Shp\\USA_States.shp")

############## get the full list of BGC zones
setwd("C:\\Colin\\Projects\\2019_CCISS")

## parameters
grid <- "WNA2"
BGC.pred.ref.WNA <- read.csv(paste("OutputData\\BGC.pred.ref", grid, "csv", sep="."))[,1]
rcp="rcp45"
proj.year=2055
BGC.pred.WNA <- read.csv(paste("OutputData\\BGC.pred", grid, rcp, proj.year, "csv", sep="."))[,1]


zone.pred.ref.WNA <- gsub("[:a-z:]","",BGC.pred.ref.WNA) 
zone.pred.ref.WNA <- gsub("[:1-9:]","",zone.pred.ref.WNA) 

zone.pred.WNA <- gsub("[:a-z:]","",BGC.pred.WNA) 
zone.pred.WNA <- gsub("[:1-9:]","",zone.pred.WNA) 

zones <- sort(unique(c(zone.pred.ref.WNA, zone.pred.WNA)))
zones <- factor(zones, levels=zones)

#BGC zones
BGCcolors <- read.csv("C:\\Colin\\Projects\\2019_CCISS\\InputData\\BGCzone_Colorscheme.csv")
colors = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)][-1]

ColScheme.WNA <- rep(NA, length(zones))
ColScheme.WNA[which(zones%in%BGCcolors$zone)] <- as.character(BGCcolors$HEX[match(zones[which(zones%in%BGCcolors$zone)], BGCcolors$zone)])
set.seed(8)
ColScheme.WNA[-which(zones%in%BGCcolors$zone)] <- sample(colors, length(zones[-which(zones%in%BGCcolors$zone)]))
ColScheme.WNA <- factor(ColScheme.WNA, levels=ColScheme.WNA)


# x11(width=6.5, height=7, pointsize=10)
png(filename=paste("Results\\CCISS.manu.BGCprojections",proj.year.focal,"png",sep="."), type="cairo", units="in", width=6.5, height=7, pointsize=10, res=600)
# mat <- matrix(c(1,1,2,3),2, byrow=T)   #define the plotting order
# layout(mat, widths=c(1,1), heights=c(1,1))   #set up the multipanel plot
par(mar=c(0.1,0.1,0.1,0.1), mgp=c(2,0.25,0))

#=============================
## Base plot
par(plt=c(0,1,0,1), bg=NA)
plot(0, col="white", xaxt="n", yaxt="n", xlab="", ylab="")
box(col="white")

#=============================
## featured GCM projection
par(plt = c(0.1,1,0.55,1), new = TRUE)
#predicted BGC zones (projected period)
rcp<- rcp.focal
proj.year <- proj.year.focal
GCM <- GCMs[4]

BGC.pred <- get(paste("BGC.pred", GCM, rcp, proj.year, sep="."))
zone.pred <- rep(NA, length(BGC.pred))
for(i in zones){ zone.pred[grep(i,BGC.pred)] <- i }
exotic <- table(zone.pred[-which(zone.pred%in%BGCcolors$zone)])
exotic <- exotic[exotic>20]
exotic <- exotic[rev(order(exotic))]
as.numeric(formatC(signif(exotic/length(zone.pred)*100,digits=3), digits=3,format="fg", flag="#"))
pred <- zone.pred
values(X) <- NA
values(X) <- factor(pred, levels=zones)[plotOrder]
values(X)[1:length(zones)] <- 1:length(zones) # this is a patch that is necessary to get the color scheme right.

plot(X, xaxt="n", yaxt="n", col=alpha(ColScheme.WNA, 1), legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE)
values(X)[-(1:length(zones))] <- NA # cover up the color bar
image(X, add=T, col="white") # cover up the color bar
values(X) <- factor(pred, levels=zones)[plotOrder] # restore the raster values
plot(bdy.bc, add=T, lwd=0.4)
# mtext(paste("(", letters[j], ") ", rcp, " ",proj.year, sep=""), side=1, line=-4.5, adj=0.12, cex=1, font=2)
mtext(paste("(B) ", GCM, "\n      (", proj.year.name[which(proj.years==proj.year)], ", " , rcp.name[which(rcps==rcp)], ")", sep=""), side=3, line=-1.5, adj=0, cex=0.8, font=2)

exotic.pct <- round(as.numeric(formatC(signif(exotic/length(zone.pred)*100,digits=3), digits=3,format="fg", flag="#")),2)
# legend(250000, 1500000, cex=0.8, title="Exotic zones", legend=paste(names(exotic), " (", exotic.pct, "%)", sep=""), fill=alpha(ColScheme.zone[as.numeric(factor(names(exotic), zone))], 1), bty="n")
legend("bottomright", cex=0.8, title="Non-BC zones\n(% of BC)", legend=paste(names(exotic), " (", exotic.pct, "%)", sep=""), fill=alpha(ColScheme.WNA[as.numeric(factor(names(exotic), zones))], 1), bty="n")

bgcs <- names(exotic)
for(bgc in bgcs){
  pts <- which(zones[values(X)]==bgc)
q=0.5
      pt <- xyFromCell(X, pts[min(which(pts >= quantile(pts, q)))])
    points(pt, pch=1, cex=1, lwd=0.8)
    text(pt-c(0, 0), bgc, pos=c(2,4,4,4,1,2)[which(bgcs==bgc)], cex=0.7, font=1, offset=0.3)
    # print(q)
  }


#=============================
## Mapped BGCs zones
par(plt = c(0, 0.45, 0.55, 0.9), new = TRUE)
#Mapped BGC zones
zones.bc <- c("BG", "BWBS", "CDF", "CWH", "ESSF", "ICH", "IDF", "MH", "MS", "PP", "SBPS", "SBS", "SWB" )
ColScheme <- as.character(BGCcolors[match(zones.bc,BGCcolors[,1]),5])

zone.bc <- rep(NA, length(BGC))
for(i in BGCcolors$zone){ zone.bc[grep(i,BGC)] <- i }
table(zone.bc)

pred <- zone.bc
pred[grep("IMA|CMA|BAFA", pred)] <- NA
X <- dem
values(X) <- NA
values(X) <- as.numeric(factor(pred))[plotOrder]
plot(X, xaxt="n", yaxt="n", col=alpha(ColScheme, 1), legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE) 
plot(bdy.bc, add=T, lwd=0.4)
par(xpd=T) 
mtext("(A) BGC zone map (1961-1990)", side=3, line=0, adj=0.05, cex=0.8, font=2)
legend("bottomleft", legend=c(zones.bc, "Alpine"), fill=c(ColScheme, "white"), bty="n", cex=0.8, inset=-0.02)
# box()


#=============================
## exotic units
par(plt = c(0.65, 1, 0.7, 1), new = TRUE)
exotic <- get(paste("BGC.pred.exotic", rcp, proj.year, sep="."))

values(X) <- NA
values(X) <- exotic[plotOrder]
breakpoints <- seq(0,1,1/length(GCMs)); length(breakpoints)
# ColScheme <- c("white", colorRampPalette(brewer.pal(9,"YlOrBr"))(length(breakpoints)-2)); length(ColScheme)
# ColScheme <- c("white", colorRampPalette(brewer.pal(9,"YlOrBr"))(length(breakpoints)-7), rep(rev(brewer.pal(9,"YlOrBr"))[1], 5)); length(ColScheme)
ColScheme <- c("white", colorRampPalette(brewer.pal(9,"Greys"))(length(breakpoints)-4), rep(rev(brewer.pal(9,"Greys"))[1], 2)); length(ColScheme)

plot(X, xaxt="n", yaxt="n", col=alpha(ColScheme, 1), legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE) 
# box(col="white", lwd=1.5)
plot(bdy.bc, add=T, lwd=0.4)
xl <- 1600000; yb <- 1000000; xr <- 1670000; yt <- 1700000
rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme, border = NA)
rect(xl,  yb,  xr,  yt)
text(rep(xr-40000,length(breakpoints))[seq(1, length(breakpoints), 3)],seq(yb,yt,(yt-yb)/(length(breakpoints)-1))[seq(1, length(breakpoints), 3)],paste(round(breakpoints*100),"%", sep="")[seq(1, length(breakpoints), 3)],pos=4,cex=0.7,font=1)
text(xl-40000, mean(c(yb,yt))-30000, paste("% Non-BC climates\nin GCM ensemble"), srt=90, pos=3, cex=0.8, font=2)
# mtext(paste(rcp, proj.year, sep="\n"), side=1, line=-3.5, adj=0.05, cex=1, font=2)
mtext(paste("(C) ", sep=""), side=3, line=-3, adj=0.05, cex=0.8, font=2)


#=============================
## single GCMs
x1 <- c(0,0.225,0,0.225)
x2 <- c(0.275,0.5,0.275,0.5)
y1 <- c(0.265, 0.265,0,0)
y2 <- c(0.505,0.505,0.24,0.24)

select <- c(8,7)
for(i in 1:2){
  GCM=GCMs[select[i]]
  rcp=rcps[1]
  proj.year=proj.year.focal
  par(plt = as.vector(rbind(x1, x2, y1, y2)[,i]), new = TRUE)
  
  BGC.pred <- get(paste("BGC.pred", GCM, rcp, proj.year, sep="."))
  zone.pred <- rep(NA, length(BGC.pred))
  for(j in zones){ zone.pred[grep(j,BGC.pred)] <- j }
  # exotic <- table(zone.pred[-which(zone.pred%in%BGCcolors$zone)])
  # exotic <- exotic[exotic>20]
  # exotic <- exotic[rev(order(exotic))]
  # as.numeric(formatC(signif(exotic/length(zone.pred)*100,digits=3), digits=3,format="fg", flag="#"))
  pred <- zone.pred
  values(X) <- NA
  values(X) <- factor(pred, levels=zones)[plotOrder]
  values(X)[1:length(zones)] <- 1:length(zones) # this is a patch that is necessary to get the color scheme right.
  
  plot(X, xaxt="n", yaxt="n", col=alpha(ColScheme.WNA, 1), legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE)
  values(X)[-(1:length(zones))] <- NA # cover up the color bar
  image(X, add=T, col="white")
  plot(bdy.bc, add=T, lwd=0.4)
  mtext(paste("(", LETTERS[4:5][i], ") ", GCM, "\n      (", proj.year.name[which(proj.years==proj.year)], ", " , c("RCP.4.5", "RCP8.5")[which(rcps==rcp)], ")", sep=""), side=3, line=-0.50, adj=0.05, cex=0.8, font=2)
  
  # exotic.pct <- round(as.numeric(formatC(signif(exotic/length(zone.pred)*100,digits=3), digits=3,format="fg", flag="#")),2)
  # legend("topright", cex=0.8, title="Exotic zones", legend=paste(names(exotic), " (", exotic.pct, "%)", sep=""), fill=alpha(ColScheme.WNA[as.numeric(factor(names(exotic), zone))], 1), bty="n")
  
  
  print(GCM)
}

#=============================
## historic
hist.year <- 2009
par(plt = as.vector(rbind(x1, x2, y1, y2)[,3]), new = TRUE)
BGC.pred <- get(paste("BGC.pred", hist.year, sep="."))
zone.pred <- rep(NA, length(BGC.pred))
for(i in zones){ zone.pred[grep(i,BGC.pred)] <- i }
# exotic <- table(zone.pred[-which(zone.pred%in%BGCcolors$zone)])
# exotic <- exotic[exotic>20]
# exotic <- exotic[rev(order(exotic))]
# as.numeric(formatC(signif(exotic/length(zone.pred)*100,digits=3), digits=3,format="fg", flag="#"))
pred <- zone.pred
values(X) <- NA
values(X) <- factor(pred, levels=zones)[plotOrder]
values(X)[1:length(zones)] <- 1:length(zones) # this is a patch that is necessary to get the color scheme right.

plot(X, xaxt="n", yaxt="n", col=alpha(ColScheme.WNA, 1), legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE)
values(X)[-(1:length(zones))] <- NA # cover up the color bar
image(X, add=T, col="white")
plot(bdy.bc, add=T, lwd=0.4)
mtext(paste("(", LETTERS[6], ") Recent (", hist.year.name[which(hist.years==hist.year)],")", sep=""), side=3, line=-0.50, adj=0.05, cex=0.8, font=2)

# exotic.pct <- round(as.numeric(formatC(signif(exotic/length(zone.pred)*100,digits=3), digits=3,format="fg", flag="#")),2)
# legend("topright", cex=0.8, title="Exotic zones", legend=paste(names(exotic), " (", exotic.pct, "%)", sep=""), fill=alpha(ColScheme.WNA[as.numeric(factor(names(exotic), zone))], 1), bty="n")

#=============================
## RCP8.5
  GCM=GCMs[4]
  rcp=rcps[2]
  proj.year=proj.years[3]
  par(plt = as.vector(rbind(x1, x2, y1, y2)[,4]), new = TRUE)
  
  BGC.pred <- get(paste("BGC.pred", GCM, rcp, proj.year, sep="."))
  zone.pred <- rep(NA, length(BGC.pred))
  for(j in zones){ zone.pred[grep(j,BGC.pred)] <- j }
  # exotic <- table(zone.pred[-which(zone.pred%in%BGCcolors$zone)])
  # exotic <- exotic[exotic>20]
  # exotic <- exotic[rev(order(exotic))]
  # as.numeric(formatC(signif(exotic/length(zone.pred)*100,digits=3), digits=3,format="fg", flag="#"))
  pred <- zone.pred
  values(X) <- NA
  values(X) <- factor(pred, levels=zones)[plotOrder]
  values(X)[1:length(zones)] <- 1:length(zones) # this is a patch that is necessary to get the color scheme right.
  
  plot(X, xaxt="n", yaxt="n", col=alpha(ColScheme.WNA, 1), legend=FALSE, legend.mar=0, maxpixels=ncell(X), bty="n", box=FALSE)
  values(X)[-(1:length(zones))] <- NA # cover up the color bar
  image(X, add=T, col="white")
  plot(bdy.bc, add=T, lwd=0.4)
  mtext(paste("(", LETTERS[7], ") ", GCM, "\n      (", proj.year.name[which(proj.years==proj.year)], ", " , c("RCP.4.5", "RCP8.5")[which(rcps==rcp)], ")", sep=""), side=3, line=-0.50, adj=0.05, cex=0.8, font=2)
  

  #=============================
  ## plot of displacement vs MAT change
  
# png(filename=paste("Results\\CCISS.BGCEDA.BGCchangeVsMATchange.png",sep="."), type="cairo", units="in", width=3.5, height=3.5, pointsize=9, res=300)
par(mar=c(2.25,3.25,0.1,0.1), mgp=c(1.75,0.25,0), plt = c(0.575, 0.999, 0.075, 0.525), new = TRUE)
plot(0, xlim=c(0,7.2), ylim=c(0,1.01), xaxs="i", yaxs="i", col="white", xaxt="n", yaxt="n",
     xlab=list(bquote(Projected~temperature~change~relative~to~"1961-90"~"("*degree*C*")"), cex=0.8), 
     ylab=list("Projected displacement of modeled BGC unit", cex=0.8))
axis(1, at=0:8, labels = 0:8, tck=0, cex.axis=0.8)
axis(2, at=seq(0,1,0.2), labels = paste(seq(0,1,0.2)*100, "%", sep=""), las=2, tck=0, cex.axis=0.8)
mtext(paste("(H) ", sep=""), side=3, line=-1.5, adj=0.05, cex=1, font=2)

select <- c(4,7,8)
int.y <- 0.07
start.y <- 0.03
pos.y <- seq(start.y,start.y+int.y*3, int.y)
offset.y <- 6
ColScheme=c("dodgerblue", "red")
for(rcp in rcps){
  for(proj.year in proj.years){
    temp.bgc <- rep(NA, length(GCMs))
    temp.zone <- rep(NA, length(GCMs))
    for(GCM in GCMs){
      temp.bgc[which(GCMs==GCM)] <- mean(get(paste("BGC.change.pred", GCM, rcp, proj.year, sep=".")) )
      temp.zone[which(GCMs==GCM)] <- mean(get(paste("zone.change.pred", GCM, rcp, proj.year, sep=".")) )
    }
    
    x <- get(paste("MAT.change", rcp, proj.year, sep="."))
    y1 <- temp.bgc
    y2 <- temp.zone
    
    points(x, y1, pch=16)
    points(x, y2, pch=1)
    
    if(rcp==rcp.focal){
      if(proj.year==proj.year.focal){
      points(x[select], y2[select], pch=21, cex=1.4, lwd=2, bg="gray")
      text(x[select], y2[select], paste(GCMs[select], " (", proj.year.name[which(proj.years==proj.year)], ", " , c("RCP.4.5", "RCP8.5")[which(rcps==rcp)], ")", sep=""), font=2, cex=0.5, pos=4)
    }
    }
    
    if(rcp==rcps[2]){
      if(proj.year==proj.years[3]){
        points(x[select[1]], y2[select[1]], pch=21, cex=1.4, lwd=2, bg="gray")
        text(x[select[1]], y2[select[1]]-0.02, paste(GCMs[select[1]], "\n(", proj.year.name[which(proj.years==proj.year)], ", " , c("RCP.4.5", "RCP8.5")[which(rcps==rcp)], ")", sep=""), font=2, cex=0.5, pos=1)
        # text(x[select[1]]+0.6, y2[select[1]]+0.025, paste(GCMs[select[1]], " (", proj.year.name[which(proj.years==proj.year)], ", " , c("RCP.4.5", "RCP8.5")[which(rcps==rcp)], ")", sep=""), font=2, cex=0.5, pos=2)
      }
    }
    
    boxplot(x, add=T, horizontal=TRUE, axes=FALSE, range=0, boxwex = 0.04, col=ColScheme[which(rcps==rcp)], at= if(rcp==rcps[1]) pos.y[which(proj.years==proj.year)]-int.y/offset.y else pos.y[which(proj.years==proj.year)]+int.y/offset.y)
    # text(max(x), ylim[1], paste(rcp.name[which(rcps==rcp)], ", ", proj.year.name[which(proj.years==proj.year)], sep=""), pos=4, cex=0.8)
    text(0, pos.y[which(proj.years==proj.year)], proj.year.name[which(proj.years==proj.year)], pos=4, cex=0.8)
    if(proj.year==proj.years[1]) text(max(x), if(rcp==rcps[1]) pos.y[which(proj.years==proj.year)]-int.y/offset.y*1.25 else pos.y[which(proj.years==proj.year)]+int.y/offset.y*1.25, 
                                      rcp.name[which(rcps==rcp)], pos=4, col=ColScheme[which(rcps==rcp)], cex=0.85, font=2)
  }
}

for(hist.year in hist.years[c(2,4)]){
  x <- get(paste("MAT.change", hist.year, sep="."))
  y1 <- mean(get(paste("BGC.change.pred", hist.year, sep=".")) )
  y2 <- mean(get(paste("zone.change.pred", hist.year, sep=".")) )
  points(x,y1, cex=1.7, pch=17)
  # text(x,y1, c("1991-\n2000", "2001-\n2010","2011-\n2017", "2017")[which(hist.years==hist.year)], pos=c(2,2,2,2)[which(hist.years==hist.year)], cex=0.7, font=2)
  points(x,y2, cex=1.4, pch=2)
  text(x,y2, hist.year.name[which(hist.years==hist.year)], pos=4, cex=0.7, font=2)
}

if(proj.year.focal==proj.years[1]){ legend(3.2, 0.48, cex=0.8, legend=c("BGC unit displacement", "BGC zone displacement", "Observed climates", ""), y.intersp = 1, pch=c(16,1,2, 1), col=c(1,1,1,"white"), pt.cex=c(1.2,1.2,1.5, 1), bty="n")
} else legend("right", cex=0.8, legend=c("BGC unit displacement", "BGC zone displacement", "Observed climates", ""), y.intersp = 1, pch=c(16,1,2, 1), col=c(1,1,1,"white"), pt.cex=c(1.2,1.2,1.5, 1), bty="n")


dev.off()



