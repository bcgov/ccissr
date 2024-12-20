
#---------------------------
#---------------------------
# testing the novelty metric for biogeoclimatic projections
#---------------------------
#---------------------------

library(climr)
library(terra)
library(data.table)
library(bcmaps)
library(ccissr)
library(ranger)
library(scales)
library(EnvStats)
library(plotly)

#---------------------------
# Data
#---------------------------

#BGC model and variable list
load("//objectstore2.nrs.bcgov/ffec/BGC_models/BGC_RFresp.Rdata") ##load RF model
pred_vars <- BGC_RFresp[["forest"]][["independent.variable.names"]] ##required predictors

# bc boundary
bc <- vect(bc_bound())
bc <- project(bc, "EPSG:4326")

# DEM
dir <- paste("//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/PRISM_dem/", sep="")
dem <- rast(paste(dir, "PRISM_dem.asc", sep=""))
dem <- aggregate(dem, fact=3)
dem <- mask(dem, bc)
dem <- trim(dem)
plot(dem)
# climate data for the biogeoclimatic projections
grid <- as.data.frame(dem, cells = TRUE, xy = TRUE)
colnames(grid) <- c("id", "lon", "lat", "elev") # rename column names to what climr expects
clim.grid <- downscale(xyz = grid,
                       gcms = list_gcms()[1],
                       ssps = list_ssps()[2],
                       gcm_periods = list_gcm_periods(),
                       run_nm = list_runs_ssp(list_gcms()[1], list_ssps()[2])[3],
                       vars = list_vars()
)
addVars(clim.grid)
clim.grid <- clim.grid[is.finite(CMD.total)] #remove NA rows to have complete cases for RF model

#historical climate for training points
pts <- fread("//objectstore2.nrs.bcgov/ffec/BGC_models/WNA_v13_50-200filtpts_15Nov.csv")
colnames(pts) <- c("id", "BGC", "lon", "lat", "elev") # rename column names to what climr expects
clim.pts <- downscale(xyz = pts,
                      vars = list_vars())
addVars(clim.pts)

# Calculate the centroid climate for the training points
clim.pts.mean <- clim.pts[, lapply(.SD, mean), by = pts$BGC, .SDcols = -c(1,2)]

# historical interannual climatic variability at the geographic centroids of the training points
pts.mean <- pts[, lapply(.SD, mean), by = BGC]
pts.mean$id <- 1:dim(pts.mean)[1]
clim.icv.pts <- downscale(xyz = pts.mean,
                          obs_years = 1951:1990,
                          obs_ts_dataset = "cru.gpcc",
                          return_refperiod = FALSE,
                          vars = list_vars())
addVars(clim.icv.pts)


#---------------------------
# EDA with Scree and 3D plots
#---------------------------

# climate data and BGC projections
clim.targets <- clim.grid[PERIOD == list_gcm_periods()[3], ]
bgc.pred <- predict(BGC_RFresp, data = clim.targets)[['predictions']]

# plots of focal analogs
pcnum <- 3
bgc.focal = "CWHxm_WA" # moderate to high novelty
bgc.focal <- "CRFdh_CA"
bgc.focal <- "ESSFwm3" # good example of true analog plus novel fringe
bgc.focal <- "SBAPcp" # example where pooled pca is necessary.
bgc.focal <- "MGPmg"
bgc.focal <- "BWBScmE" # good example of novelty along a low-variance dimension.
bgc.focal <- "IDFmw2" # significant separation in the 4th PC.
bgc.focal <- "CWHvh1" # shows an apparent bug in the ICV centroid.
bgc.focal <- "IDFdxx_WY" #
bgc.focal <- "IDFdk5" # low novelty in pred_vars but high novelty in basic variables
analog_novelty(clim.targets = clim.targets,
               clim.analogs = clim.pts,
               label.targets = bgc.pred,
               label.analogs = pts$BGC,
               # vars = pred_vars[-which(pred_vars=="CMI")], # remove CMI as it is NA along the coast (climr bug)
               vars = as.vector(outer(c("Tmin", "Tmax", "PPT"), c("wt", "sp", "sm", "at"), paste, sep = "_")),
               pcs = 3,
               analog.focal = bgc.focal,
               # logVars = FALSE,
               plotScree = TRUE,
               clim.icvs <- clim.icv.pts,
               label.icvs <- pts.mean$BGC[clim.icv.pts$id],
               # plot2d = TRUE,
               plot3d = TRUE,
               plot3d.pcs=c(1,2,3), 
               plot3d.candidates = TRUE
)


pca <- prcomp(clim.pts, scale=TRUE)


#---------------------------
# Novelty maps
#---------------------------

# Color Scheme for sigma novelty
breakseq <- c(0,4,8)
breakpoints <- c(seq(breakseq[1], breakseq[3], 0.01),199); length(breakpoints)
ColScheme <- c(colorRampPalette(c("gray90", "gray50", "#FFF200", "#CD0000", "black"))(length(breakpoints)))

# template Raster
X <- dem
values(X) <- NA

# climate data and BGC projections
clim.targets <- clim.grid[PERIOD == list_gcm_periods()[1], ]
bgc.pred <- predict(BGC_RFresp, data = clim.targets)[['predictions']]

par(mar=c(1,1,1,1), mfrow=c(1,1))
novelty <- analog_novelty(clim.targets = clim.targets,
                          clim.analogs = clim.pts,
                          label.targets = bgc.pred,
                          label.analogs = pts$BGC,
                          # vars = pred_vars[-which(pred_vars=="CMI")], # remove CMI as it is NA along the coast (climr bug)
                          vars = as.vector(outer(c("Tmin", "Tmax", "PPT"), c("wt", "sp", "sm", "at"), paste, sep = "_")),
                          clim.icvs <- clim.icv.pts,
                          label.icvs <- pts.mean$BGC[clim.icv.pts$id],
                          weight.icv = 0.5,
                          threshold = 0.95,
                          pcs = NULL, 
                          logVars = TRUE
                          
)
X[clim.targets[, id]] <- novelty
plot(X, col=ColScheme, axes=F)
mtext("Sigma novelty", side=4, line=-4.5, adj = 0.5, font=2)

# Maps at different ICV weights
par(mar=c(1,1,1,1), mfrow=c(2,2))
weight.icvs=c(0, 0.33, 0.66, 1)
for(weight.icv in weight.icvs){
  novelty <- analog_novelty(clim.targets = clim.targets,
                            clim.analogs = clim.pts,
                            label.targets <- bgc.pred,
                            label.analogs <- pts$BGC,
                            clim.icvs <- clim.icv.pts,
                            label.icvs <- pts.mean$BGC[clim.icv.pts$id],
                            weight.icv = weight.icv,
                            # vars = pred_vars[-which(pred_vars=="CMI")],
                            vars = as.vector(outer(c("Tmin", "Tmax", "PPT"), c("wt", "sp", "sm", "at"), paste, sep = "_")),
                            threshold = 0.95,
                            pcs = NULL
  )
  X[clim.targets[, id]] <- novelty
  plot(X, col=ColScheme, axes=F)
  mtext("Sigma novelty", side=4, line=-4.5, adj = 0.5, font=2)
  mtext(paste0("RF variables", "\n", "ICV weight = ", weight.icv), line=-3.5, adj = 0.975, )
}

# Maps at different ICV weights
par(mar=c(1,1,1,1), mfrow=c(2,2))
weight.icvs=c(0, 0.33, 0.66, 1)
for(weight.icv in weight.icvs){
  novelty <- analog_novelty(clim.targets = clim.targets,
                            clim.analogs = clim.pts,
                            label.targets <- bgc.pred,
                            label.analogs <- pts$BGC,
                            clim.icvs <- clim.icv.pts,
                            label.icvs <- pts.mean$BGC[clim.icv.pts$id],
                            weight.icv = weight.icv,
                            # vars = pred_vars[-which(pred_vars=="CMI")],
                            vars = as.vector(outer(c("Tmin", "Tmax", "PPT"), c("wt", "sp", "sm", "at"), paste, sep = "_")),
                            threshold = 0.95,
                            pcs = NULL
  )
  X[clim.targets[, id]] <- novelty
  plot(X, col=ColScheme, axes=F)
  mtext("Sigma novelty", side=4, line=-4.5, adj = 0.5, font=2)
  mtext(paste0("RF variables", "\n", "ICV weight = ", weight.icv), line=-3.5, adj = 0.975, )
}

# Map at increasing PCs
par(mar=c(1,1,1,1), mfrow=c(2,2))
for(pcnum in c(2,3,4,5)){
  novelty <- analog_novelty(clim.targets = clim.targets,
                            clim.analogs = clim.pts,
                            label.targets = bgc.pred,
                            label.analogs = pts$BGC,
                            # vars = pred_vars[-which(pred_vars=="CMI")],
                            vars = as.vector(outer(c("Tmin", "Tmax", "PPT"), c("wt", "sp", "sm", "at"), paste, sep = "_")),
                            clim.icvs <- clim.icv.pts,
                            label.icvs <- pts.mean$BGC[clim.icv.pts$id],
                            weight.icv = 0.5,
                            threshold = 0.95,
                            pcs = pcnum
                            
  )
  X[clim.targets[, id]] <- novelty
  plot(X, col=ColScheme, axes=F)
  mtext("Sigma novelty", side=4, line=-4.5, adj = 0.5, font=2)
  mtext(paste0("RF variables", "\n", pcnum, "PCs"), line=-3.5, adj = 0.975, )
}

# Maps with and without log-transformation
clim.targets <- clim.grid[PERIOD == list_gcm_periods()[1], ]
par(mar=c(1,1,1,1), mfrow=c(1,2))
for(log in c(TRUE, FALSE)){
  novelty <- analog_novelty(clim.targets = clim.targets,
                            clim.analogs = clim.pts,
                            label.targets <- bgc.pred,
                            label.analogs <- pts$BGC,
                            clim.icvs <- clim.icv.pts,
                            label.icvs <- pts.mean$BGC[clim.icv.pts$id],
                            logVars = log,
                            # vars = pred_vars[-which(pred_vars=="CMI")],
                            vars = as.vector(outer(c("Tmin", "Tmax", "PPT"), c("wt", "sp", "sm", "at"), paste, sep = "_")),
                            threshold = 0.95,
                            pcs = NULL
  )
  X[clim.targets[, id]] <- novelty
  plot(X, col=ColScheme, axes=F)
  mtext("Sigma novelty", side=4, line=-4.5, adj = 0.5, font=2)
  mtext(paste0("RF variables", "\n", "logVars = ", log), line=-3.5, adj = 0.975, )
}

#---------------------------
# novelty over time
#---------------------------

par(mar=c(1,1,1,1), mfrow=c(2,2))
for(i in 0:3){ # i is the iteration through gcm_periods
  
  # subset climate data to selected time period
  if(i==0){
    clim.targets <- clim.grid[PERIOD == "1961_1990", ]
  } else clim.targets <- clim.grid[PERIOD == list_gcm_periods()[i], ]
  
  # BGC projections
  bgc.pred <- predict(BGC_RFresp, data = clim.targets)[['predictions']]
  
  novelty <- analog_novelty(clim.targets = clim.targets,
                            clim.analogs = clim.pts,
                            label.targets = bgc.pred,
                            label.analogs = pts$BGC,
                            vars = pred_vars[-which(pred_vars=="CMI")],
                            clim.icvs <- clim.icv.pts,
                            label.icvs <- pts.mean$BGC[clim.icv.pts$id]
  )
  X[clim.targets[, id]] <- novelty
  plot(X, col=ColScheme, axes=F)
  mtext("Sigma novelty", side=4, line=-4.5, adj = 0.5, font=2)
  mtext(paste0("RF variables", "\n", list_gcm_periods()[i]), line=-3.5, adj = 0.975, )
}

dat <- rast("../../Novelty_EC-Earth3_2061_2080.tif")
