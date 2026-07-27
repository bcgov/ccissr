
## -------------------------------------------------
## -------------------------------------------------
## bivariate kernel density distributions of baseline and projected bgc projections
## -------------------------------------------------
## -------------------------------------------------


library(climr)
library(ccissr)
library(data.table)
library(terra)
library(leaflet)
library(ranger) # For RF
library(MASS)
library(scales)

## -------------------------------------------------
## bgc metadata
data("WNA_BGCs")


## -------------------------------------------------
## Input data (need to use 250m DEM because otherwise error in elevation values is too high)

# dir <- "//objectstore2.nrs.bcgov/ffec/BGC_models/" 
dir <- "C:/Users/CMAHONY/Data/BGC_models/" #local copy, for speed

BGCmodel <- readRDS(paste0(dir, "BGCmodel_WNA_V4.2gini.rds")) #Kiri trained this on Thufir using the Gini split rule. 

studyname <- "BC"
bdy <- vect("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/bdy.BC.shp")
bdy <- project(bdy, "EPSG:4326")
dem <- rast("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling//WNA_DEM_4326_clipped.tif")
dem <- crop(dem, ext(bdy) )
# dem <- aggregate(dem, fact=5)
dem <- mask(dem, bdy) 
# plot(dem)
X <- dem # template raster for testing
values(X) <- NA

bgcs <- rast("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/WNA_BGC_v13_13Apr2026_raster.tif")
bgcs <- crop(bgcs, ext(bdy))
bgcs <- mask(bgcs, dem) 
bgc_levels <- cats(bgcs)[[1]]  # Extract the category mapping
vals_char <- bgc_levels$BGC[match(values(bgcs), bgc_levels$value)] # Convert numeric values to character labels
values(bgcs) <- factor(vals_char, levels = WNA_BGCs$BGC) # refactor with full set of levels

## study area points
points <- as.data.table(dem, cells=T, xy=T)
colnames(points) <- c("id", "lon", "lat", "elev")
points <- points[,c(2,3,4,1)] #restructure for climr input
values(X) <- NA; values(X)[points$id] <- points$el ; plot(X)

## add bgc label to points
points.bgcs <- as.data.table(bgcs, cells=T, xy=T)
colnames(points.bgcs) <- c("id", "lon", "lat", "bgc")
points[points.bgcs[, .(id, bgc)], bgc := i.bgc, on = "id"]
points[, bgc.zone := WNA_BGCs[match(points$bgc, WNA_BGCs$BGC), Zone]]
values(X) <- NA; values(X)[points$id] <- points$zone ; plot(X)

#take a sample of 20000 points from each bgc zone to reduce the size of the file (for computational feasibility)
points.sampled <- points[, .SD[sample(.N, min(.N, 20000))], by = zone] 
points.sampled[,id := 1:(dim(points.sampled)[1])]
dim(points.sampled)

## -------------------------------------------------
## climate data for all points

clim_ref <- downscale(
  xyz = points.sampled,
  which_refmap = "refmap_climr",
  return_refperiod = TRUE, # Also return the 1961-1990 normals period.
  vars = list_vars()
)
ccissr::addVars(clim_ref)
clim_ref[!is.finite(CMD.total), CMD.total := 0]
clim_ref[!is.finite(Eref_sp), Eref_sp := 0]

clim_proj <- downscale(
  xyz = points.sampled,
  gcms = list_gcms()[5],
  ssps = list_ssps()[2],
  gcm_periods = list_gcm_periods()[3],
  run_nm = list_runs_ssp(list_gcms()[5], list_ssps()[2])[4],
  which_refmap = "refmap_climr",
  return_refperiod = FALSE, # Also return the 1961-1990 normals period.
  vars = list_vars()
)
ccissr::addVars(clim_proj)
clim_proj[!is.finite(CMD.total), CMD.total := 0]
clim_proj[!is.finite(Eref_sp), Eref_sp := 0]

## -------------------------------------------------
## predictions

## V4.2gini model predictions for reference period
preds_ref_vec <- predict(BGCmodel, data = clim_ref)$prediction

## V4.2gini model predictions for future period
preds_proj_vec <- predict(BGCmodel, data = clim_proj)$prediction


## -------------------------------------------------
## bivariate kernel density distributions of baseline and future latitude and elevation for each BGC zone
## -------------------------------------------------

zone.ref <- WNA_BGCs[match(preds_ref_vec, WNA_BGCs$BGC), Zone]

zone.proj <- WNA_BGCs[match(preds_proj_vec, WNA_BGCs$BGC), Zone]

zones <- unique(WNA_BGCs[DataSet=="BC", Zone ])

dir <- "C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/CCISS_Manual/Figures/"
png(filename=paste0(dir, "/ElevShift.zones.png"), type="cairo", units="in", width=6.5, height=6.5, pointsize=10, res=300)

mat <- matrix(c(17, 1:4, 17, 5:8, 17, 9:12, 17, 13:16, rep(18,5)),5, byrow=T)   #define the plotting order
layout(mat, widths=c(0.1,1,1,1,1), heights=c(1,1,1,1,0.1))   #set up the multipanel plot
par(mar=c(1.5,1.5,0.8,0.8), mgp=c(1.5,0.25,0), tck= -0.01)

for(zone in zones){
  elev.ref <- points.sampled[which(zone.ref==zone), elev]
  lat.ref <- points.sampled[which(zone.ref==zone), lat]
  elev.proj <- points.sampled[which(zone.proj==zone), elev]
  lat.proj <- points.sampled[which(zone.proj==zone), lat]
  
  x <- c(elev.ref, elev.proj)
  y <- c(lat.ref, lat.proj)
  
  # KDE parameters
  h.factor <- 0.1
  prob.threshold <- 0.95 #probability contour to draw as polygon
  
  # Kernel density estimation for ref period
  k <- kde2d(elev.ref, lat.ref, n=500, h=c(h.factor*diff(range(x)), h.factor*diff(range(y))), lims = c(range(x)+c(-200, 200), range(y)+c(-2, 20)))
  dx <- diff(k$x[1:2]); dy <- diff(k$y[1:2]) # grid cell area
  z <- k$z / sum(k$z * dx * dy) # normalize density surface
  zvec <- as.vector(z); # flatten
  zsort <- sort(zvec, decreasing = TRUE) # sort
  cumprob <- cumsum(zsort * dx * dy) # cumulative probability mass
  lev <- zsort[min(which(cumprob >= prob.threshold))] # threshold for specified probability threshold HDR
  cl.ref <- contourLines(x = k$x, y = k$y, z = k$z, levels = lev)
  
  # Kernel density estimation for proj period
  k <- kde2d(elev.proj, lat.proj, n=500, h=c(h.factor*diff(range(x)), h.factor*diff(range(y))), lims = c(range(x)+c(-200, 200), range(y)+c(-2, 20)))
  dx <- diff(k$x[1:2]); dy <- diff(k$y[1:2]) # grid cell area
  z <- k$z / sum(k$z * dx * dy) # normalize density surface
  zvec <- as.vector(z); # flatten
  zsort <- sort(zvec, decreasing = TRUE) # sort
  cumprob <- cumsum(zsort * dx * dy) # cumulative probability mass
  lev <- zsort[min(which(cumprob >= prob.threshold))] # threshold for specified probability threshold HDR
  cl.proj <- contourLines(x = k$x, y = k$y, z = k$z, levels = lev)
  
  # all x coordinates
  xrange <- range(c(
    unlist(lapply(cl.ref, `[[`, "x")),
    unlist(lapply(cl.proj, `[[`, "x"))
  ))
  
  # all y coordinates
  yrange <- range(c(
    unlist(lapply(cl.ref, `[[`, "y")),
    unlist(lapply(cl.proj, `[[`, "y"))
  ))
  
  # plot(x,y, xaxs="i", yaxs="i", col="white")
  plot(1, type = "n", xlab = "", ylab="", xlim = range(xrange), ylim = range(yrange))
  
  # plot polygons
  for(i in 1:length(cl.ref)){polygon(cl.ref[[i]]$x, cl.ref[[i]]$y, border = "grey", lwd = 1, col = adjustcolor("grey", alpha.f = 0.2))}   # draw polygon
  for(i in 1:length(cl.proj)){polygon(cl.proj[[i]]$x, cl.proj[[i]]$y, border = "dodgerblue", lwd = 1, col = adjustcolor("dodgerblue", alpha.f = 0.2))}   # draw polygon
  # s <- sample(1:length(elev.ref), 10000)
  # points.sampled(elev.ref,lat.ref, col=alpha("dodgerblue", 0.25), pch=16, cex=0.5)
  
  mtext(paste0("(", letters[which(zones==zone)], ") ", zone), line=-1.5, side=3, adj=0.025, cex=0.8, font=1)
  
  if(zone==zones[1]){
    legend("bottomleft", legend=c("Baseline (1961-1990)", "Future (2041-2060)"), fill=alpha(c("grey", "dodgerblue"), 0.2), border=c("grey", "dodgerblue"), bty="n")
  }
  
  # plot(x,y, col="white", xlab = "Elevation (m)", ylab="Latitude (deg.)")
  # points.sampled(elev.ref,lat.ref, col=alpha("dodgerblue", 0.25), pch=16, cex=1.5)
  # points.sampled(elev.proj,lat.proj, col=alpha("grey40", 0.85), cex=1.5)
  
  print(zone)
  
}
par(mar=c(0,0,0,0))

plot(1, type="n", axes=F, xlab="", ylab="")  
text(1,1, "Latitude (degrees)", srt=90, font=1,cex=1.3)  

par(mar=c(0,0,0,0))
plot(1, type="n", axes=F, xlab="", ylab="")  
text(1,1, "Elevation (m)", font=1,cex=1.3)  

dev.off()




## -------------------------------------------------
## bivariate kernel density distributions of bgc map and baseline prediction latitude and elevation for each BGC zone
## -------------------------------------------------

zone.ref <- WNA_BGCs[match(preds_ref_vec, WNA_BGCs$BGC), Zone]

zone.map <- points.sampled$bgc.zone

zones <- unique(WNA_BGCs[DataSet=="BC", Zone ])

dir <- "C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/CCISS_Manual/Figures/"
png(filename=paste0(dir, "/ElevShift.zones.baseline.png"), type="cairo", units="in", width=6.5, height=6.5, pointsize=10, res=300)

mat <- matrix(c(17, 1:4, 17, 5:8, 17, 9:12, 17, 13:16, rep(18,5)),5, byrow=T)   #define the plotting order
layout(mat, widths=c(0.1,1,1,1,1), heights=c(1,1,1,1,0.1))   #set up the multipanel plot
par(mar=c(1.5,1.5,0.8,0.8), mgp=c(1.5,0.25,0), tck= -0.01)

for(zone in zones){
  elev.ref <- points.sampled[which(zone.ref==zone), elev]
  lat.ref <- points.sampled[which(zone.ref==zone), lat]
  elev.map <- points.sampled[which(zone.map==zone), elev]
  lat.map <- points.sampled[which(zone.map==zone), lat]
  
  x <- c(elev.ref, elev.map)
  y <- c(lat.ref, lat.map)
  
  # KDE parameters
  h.factor <- 0.1
  prob.threshold <- 0.95 #probability contour to draw as polygon
  
  # Kernel density estimation for ref period prediction
  k <- kde2d(elev.ref, lat.ref, n=500, h=c(h.factor*diff(range(x)), h.factor*diff(range(y))), lims = c(range(x)+c(-200, 200), range(y)+c(-2, 20)))
  dx <- diff(k$x[1:2]); dy <- diff(k$y[1:2]) # grid cell area
  z <- k$z / sum(k$z * dx * dy) # normalize density surface
  zvec <- as.vector(z); # flatten
  zsort <- sort(zvec, decreasing = TRUE) # sort
  cumprob <- cumsum(zsort * dx * dy) # cumulative probability mass
  lev <- zsort[min(which(cumprob >= prob.threshold))] # threshold for specified probability threshold HDR
  cl.ref <- contourLines(x = k$x, y = k$y, z = k$z, levels = lev)
  
  # Kernel density estimation for map 
  k <- kde2d(elev.map, lat.map, n=500, h=c(h.factor*diff(range(x)), h.factor*diff(range(y))), lims = c(range(x)+c(-200, 200), range(y)+c(-2, 20)))
  dx <- diff(k$x[1:2]); dy <- diff(k$y[1:2]) # grid cell area
  z <- k$z / sum(k$z * dx * dy) # normalize density surface
  zvec <- as.vector(z); # flatten
  zsort <- sort(zvec, decreasing = TRUE) # sort
  cumprob <- cumsum(zsort * dx * dy) # cumulative probability mass
  lev <- zsort[min(which(cumprob >= prob.threshold))] # threshold for specified probability threshold HDR
  cl.map <- contourLines(x = k$x, y = k$y, z = k$z, levels = lev)
  
  # all x coordinates
  xrange <- range(c(
    unlist(lapply(cl.ref, `[[`, "x")),
    unlist(lapply(cl.map, `[[`, "x"))
  ))
  
  # all y coordinates
  yrange <- range(c(
    unlist(lapply(cl.ref, `[[`, "y")),
    unlist(lapply(cl.map, `[[`, "y"))
  ))
  
  # plot(x,y, xaxs="i", yaxs="i", col="white")
  plot(1, type = "n", xlab = "", ylab="", xlim = range(xrange), ylim = range(yrange))
  
  # plot polygons
  for(i in 1:length(cl.ref)){polygon(cl.ref[[i]]$x, cl.ref[[i]]$y, border = "grey", lwd = 1, col = adjustcolor("grey", alpha.f = 0.2))}   # draw polygon
  for(i in 1:length(cl.map)){polygon(cl.map[[i]]$x, cl.map[[i]]$y, border = "dodgerblue", lwd = 1, col = adjustcolor("dodgerblue", alpha.f = 0.2))}   # draw polygon
  # s <- sample(1:length(elev.ref), 10000)
  # points.sampled(elev.ref,lat.ref, col=alpha("dodgerblue", 0.25), pch=16, cex=0.5)
  
  mtext(paste0("(", letters[which(zones==zone)], ") ", zone), line=-1.5, side=3, adj=0.025, cex=0.8, font=1)
  
  if(zone==zones[1]){
    legend("bottomleft", legend=c("Baseline prediction", "BGC map"), fill=alpha(c("grey", "dodgerblue"), 0.2), border=c("grey", "dodgerblue"), bty="n")
  }
  
  # plot(x,y, col="white", xlab = "Elevation (m)", ylab="Latitude (deg.)")
  # points.sampled(elev.ref,lat.ref, col=alpha("dodgerblue", 0.25), pch=16, cex=1.5)
  # points.sampled(elev.map,lat.map, col=alpha("grey40", 0.85), cex=1.5)
  
  print(zone)
  
}
par(mar=c(0,0,0,0))

plot(1, type="n", axes=F, xlab="", ylab="")  
text(1,1, "Latitude (degrees)", srt=90, font=1,cex=1.3)  

par(mar=c(0,0,0,0))
plot(1, type="n", axes=F, xlab="", ylab="")  
text(1,1, "Elevation (m)", font=1,cex=1.3)  

dev.off()















