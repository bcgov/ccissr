### CCISS Dashboard Plots
### Kiri Daust

devtools::load_all()
library(terra)
library(data.table)
library(RColorBrewer)
library(climr)
library(ggplot2)

dem2 <- rast("dem_BC2kmGrid.tif")
dem_table <- climr::dem_to_table(dem2)

BGCmodel <- readRDS("../Common_Files/WNA_BGCv13/BGCmodel_WNA_V4.2gini.rds")

bgcs <- vect("../Common_Files/WNA_BGCv13/BGC_BC_v13_Jan21.gdb/")
bgcs <- bgcs["MAP_LABEL"]
names(bgcs) <- "BGC"
bgc_template <- make_bgc_template(dem2, bgcs)

con <- dbCon_cciss("./bc_2km_v13_26.duckdb", threads = 8)
dbPopulate(con, bgc_template)

vars_needed <- c("CMD_sm", "DDsub0_sp", "DD5_sp", "Eref_sm", "Eref_sp", "EXT", 
                 "MWMT", "NFFD_sm", "NFFD_sp", "PAS_an", "PAS_sp", "SHM", "Tave_sm", "Tave_sp", "Tmax_sm", "Tmax_sp", "Tmin_an", "Tmin_at", "Tmin_sm", 
                 "Tmin_sp", "Tmin_wt","CMI_an", "PPT_05","PPT_06","PPT_07","PPT_08","PPT_09","PPT_at","PPT_wt","CMD_07","CMD_an"
)

gcms_cciss <- c("ACCESS-ESM1-5", "CNRM-ESM2-1", "EC-Earth3", "GFDL-ESM4",
                "GISS-E2-1-G", "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0")

predict_bgc(con, dem_table, BGCmodel, vars_needed, gcms_cciss, periods_use = list_gcm_periods(), 
            max_runs_use = 3L, obs_2001_2020 = TRUE)
summarise_preds(con)

siteseries_preds(con)

spp_list <- c("Pl","Sx","Fd") #
for(spp in spp_list){
  message("Processing ", spp)
  cciss_full_species(con, spp)
}

##species bubbleplots
perexp <- spp_persist_expand(con, spp_list = c("Pl","Sx","Fd"), fractional = FALSE)
spp_bubbleplot(perexp, period = "2081_2100", scenario = "ssp245", species.focal = "Fd", edatope = c("C4"))

##alluvial plots
plot_spparea(con, spp = "Fd", edatope = "C4", fractional = FALSE, save_png = FALSE)

## 2 panel maps
bc_ol <- vect("data-raw/data_tables/bc_outline.gpkg")
plot_SuitabilityChangeMap(con, bgc_template = bgc_template, outline = bc_ol, spp = "Fd", edatope = "C4", period = "2001_2020", three_panel = FALSE, save_png = FALSE) 

##reference suit maps
plot_reference_suitability(con, dem2, bc_ol, spp = "Fd", edatope = "C4")
