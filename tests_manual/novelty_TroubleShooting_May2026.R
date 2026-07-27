
#---------------------------
#---------------------------
# the novelty metric is producing erroneous results, notably high novelty on the coast even in the baseline predictions
# This is a troubleshooting script. 
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

# Color Scheme for sigma novelty
breakseq <- c(0,4,8)
breakpoints <- c(seq(breakseq[1], breakseq[3], 0.01),199); length(breakpoints)
ColScheme <- c(colorRampPalette(c("gray90", "gray50", "#FFF200", "#CD0000", "black"))(length(breakpoints)))


#---------------------------
# Data
#---------------------------

#BGC model and variable list
# BGCmodel <- readRDS("//objectstore2.nrs.bcgov/ffec/BGC_models/BGCmodel_WNA_V2.1.rds") # permanent storage but too slow
BGCmodel <- readRDS("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/Trained_Models//BGCmodel_WNA_V2.1.rds")
pred_vars <- BGCmodel[["forest"]][["independent.variable.names"]] ##required predictors

# bc boundary
bc <- vect(bc_bound())
bc <- project(bc, "EPSG:4326")

# DEM
# dir <- paste("//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/PRISM_dem/", sep="") # permanent storage but too slow
dir <- paste("C:/Users/CMAHONY/OneDrive - Government of BC/Data/PRISM_dem/", sep="")
dem <- rast(paste(dir, "PRISM_dem.asc", sep=""))
dem <- aggregate(dem, fact=3)
dem <- mask(dem, bc)
dem <- trim(dem)
plot(dem)

# template Raster
X <- dem
values(X) <- NA

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
setnames(clim.grid,  new = c("PAS", "Tmin", "CMI"),  old = c("PAS_an", "Tmin_an", "CMI_an"))

#historical climate for training points
pts <- fread("//objectstore2.nrs.bcgov/ffec/BGC_models/points_WNA_simple200.csv")
clim.pts <- downscale(xyz = pts,
                      vars = list_vars())
addVars(clim.pts)
clim.pts <- pts[clim.pts, on = "id"]
setnames(clim.pts,  new = c("PAS", "Tmin", "CMI"),  old = c("PAS_an", "Tmin_an", "CMI_an"))

# Calculate the centroid climate for the training points
clim.pts.mean <- clim.pts[, lapply(.SD, mean), by = BGC, .SDcols = -c("id", "PERIOD")]

# historical interannual climatic variability at the geographic centroids of the training points
pts.mean <- pts[, lapply(.SD, mean), by = BGC]
pts.mean$id <- 1:dim(pts.mean)[1]
clim.icv.pts <- downscale(xyz = pts.mean,
                          obs_years = 1951:1990,
                          obs_ts_dataset = "cru.gpcc",
                          return_refperiod = FALSE,
                          vars = list_vars())
addVars(clim.icv.pts)
clim.icv.pts <- pts.mean[clim.icv.pts, on = "id"]
setnames(clim.icv.pts,  new = c("PAS", "Tmin", "CMI"),  old = c("PAS_an", "Tmin_an", "CMI_an"))



#---------------------------
# Novelty of baseline prediction
#---------------------------


clim.targets <- clim.grid[PERIOD == "1961_1990", ]
bgc.pred <- predict(BGCmodel, data = clim.targets)[['predictions']]

novelty <- analog_novelty(clim.targets = clim.targets,
                          clim.analogs = clim.pts,
                          label.targets = bgc.pred,
                          label.analogs = pts$BGC,
                          vars = as.vector(outer(c("Tmin", "Tmax", "PPT"), c("wt", "sp", "sm", "at"), paste, sep = "_")),
                          clim.icvs <- clim.icv.pts,
                          label.icvs = clim.icv.pts$BGC
)

par(mar=c(1,1,1,1), mfrow=c(1,1))
X[clim.targets[, id]] <- novelty
plot(X, col=ColScheme, axes=F)
mtext("Sigma novelty", side=4, line=-4.5, adj = 0.5, font=2)
mtext(paste0("basic variables", "\n", "1961_1990"), line=-3.5, adj = 0.975, )


#---------------------------
# focal bgc
#---------------------------

bgc.focal = "CWHwh1" # moderate to high novelty

novelty <- analog_novelty(clim.targets = clim.targets,
                          clim.analogs = clim.pts,
                          label.targets = bgc.pred,
                          label.analogs = pts$BGC,
                          vars = as.vector(outer(c("Tmin", "Tmax", "PPT"), c("wt", "sp", "sm", "at"), paste, sep = "_")),
                          clim.icvs = clim.icv.pts,
                          label.icvs = clim.icv.pts$BGC,
                          analog.focal = bgc.focal,
                          plotScree = TRUE,
                          plot2d = TRUE,
                          plot3d = TRUE,
                          plot3d.pcs=c(1,2,3), 
                          plot3d.candidates = FALSE
                          )


# map of focal bgc predictions and analog sample
values(X) <- NA
X[clim.targets[, id]] <- 1
plot(X)
# analog sample: 
points(pts[BGC==bgc.focal,.(lon,lat)], pch=16, cex=0.6)
# predicted distribution: 
points(grid[which(bgc.pred==bgc.focal),2:3], col="blue", pch=1, cex=0.6)
# the spatial distribution is the same 

# are the elevations correct? 
x <- pts[BGC==bgc.focal,elev]
y <- grid[which(bgc.pred==bgc.focal),4]
hist(x)
hist(y)
# yes

# compare climates 
clim.pts.focal <- clim.pts[which(pts$BGC==bgc.focal),.(Tmin_sm, PPT_sm)]
clim.targets.focal <- clim.targets[which(bgc.pred==bgc.focal),.(Tmin_sm, PPT_sm)]
par(mar=c(3,3,1,1))
plot(clim.pts.focal, col="blue")
points(clim.targets.focal, pch=16)
legend("topleft", legend=c("analog", "targets"), pch=c(1,16), col=c("blue", 1))
# the analog climates are messed up. 

# try downscaling the analog climates just for the focal analog
pts.focal <- pts[BGC==bgc.focal,]
clim.pts.focal <- downscale(xyz = pts.focal, vars = list_vars())
points(clim.pts.focal[,.(Tmin_sm, PPT_sm)], pch=16, col="red")
legend("topleft", legend=c("targets", "analog - bulk downscale", "analog - focal downscale"), pch=c(16,1,16), col=c(1, "blue", "red"))
# this returns the correct climate. it seems like the rows are getting scrambled during downscale of the bulk file. 


