### CCISS Dashboard Plots
### Kiri Daust
### with VC edits July 2026

# Need to install climr version 0.2.23 from the devl branch
#remotes::install_github("bcgov/climr@devl")

# Build ccissr
devtools::load_all()

library(terra)
library(data.table)
library(RColorBrewer)
library(climr)
#packageVersion("climr") #check the package version
library(ggplot2)
library(DBI)
library(duckdb)

dem2 <- rast("C:/Vanessa/GitHub/CCISSfiles/dem_BC2kmGrid.tif")
dem_table <- climr::dem_to_table(dem2)

BGCmodel <- readRDS("C:/Vanessa/GitHub/CCISSfiles/BGCmodel_WNA_V4.2gini.rds")

## spatial file of BGCs
bgcs <- vect("C:/Vanessa/GitHub/CCISSfiles/BGC_BC_v13_Jan21.gdb/")
bgcs <- bgcs["MAP_LABEL"]
names(bgcs) <- "BGC"
bgc_template <- make_bgc_template(dem2, bgcs)

#### Kiri sent me the duckdb file
con <- dbCon_cciss("C:/Vanessa/GitHub/CCISSfiles/bc_2km_July7.duckdb", threads = 8) 
dbExecute(con, "SET memory_limit = '16GB'") ##adjust this if you run out of memory
dbExecute(con, "SET threads = 8")
#having issues where R can not find the temporary directory, create a temp directory
dbExecute(con, "SET temp_directory = 'C:/Vanessa/GitHub/CCISSfiles/temp'")
dir.create("C:/Vanessa/GitHub/CCISSfiles/temp", showWarnings = FALSE)

############################## SKIP this next part #################################

## create/open duckdb file
con <- dbCon_cciss("./bc_2km_v13_26.duckdb", threads = 8) ##name this whatever you want
dbExecute(con, "SET memory_limit = '16GB'") ##adjust this if you run out of memory
dbExecute(con, "SET threads = 8")

dbPopulate(con, bgc_template)

vars_needed <- c("CMD_sm", "DDsub0_sp", "DD5_sp", "Eref_sm", "Eref_sp", "EXT", 
                 "MWMT", "NFFD_sm", "NFFD_sp", "PAS_an", "PAS_sp", "SHM", "Tave_sm", "Tave_sp", "Tmax_sm", "Tmax_sp", "Tmin_an", "Tmin_at", "Tmin_sm", 
                 "Tmin_sp", "Tmin_wt","CMI_an", "PPT_05","PPT_06","PPT_07","PPT_08","PPT_09","PPT_at","PPT_wt","CMD_07","CMD_an"
)

gcms_cciss <- c("ACCESS-ESM1-5", "CNRM-ESM2-1", "EC-Earth3", "GFDL-ESM4",
                "GISS-E2-1-G", "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0")

## this will take a couple of hours
predict_bgc(con, dem_table, BGCmodel, vars_needed, gcms_cciss, periods_use = list_gcm_periods(), 
            max_runs_use = 3L, obs_2001_2020 = TRUE)

summarise_preds(con)

###################### continue here ##################################################


## siteseries predictions (edatopic overlap)
siteseries_preds(con)

## create cciss species predictions in database
spp_list <- c("Pl","Sx","Fd","Cw","Hw","Py", "Bl","At", "Ac", "Ep", "Yc", "Pw", "Ss", "Bg", "Lw", "Mb", "Ba")
for(spp in spp_list){
  message("Processing ", spp)
  cciss_full_species(con, spp)
}

#### Test figures for Fd C4 ####

## calculate persistence and expansion
perexp <- spp_persist_expand(con, spp_list = spp_list, fractional = FALSE) #fractional spp suit = FALSE (binary suit yes or no)

## species bubbleplots
spp_bubbleplot(perexp, period = "2041_2060", scenario = "ssp245", species.focal = "Fd", edatope = c("C4"))

## alluvial plots - can set save_png to true when you're ready
plot_spparea(con, spp = "Fd", edatope = "C4", fractional = FALSE, save_png = FALSE)

## 2 panel maps
bc_ol <- vect("data-raw/data_tables/bc_outline.gpkg")
plot_SuitabilityChangeMap(con, bgc_template = bgc_template, outline = bc_ol, spp = "Fd", edatope = "C4", period = "2021_2040", three_panel = FALSE, save_png = FALSE) 

## reference suit maps
plot_reference_suitability(con, dem2, bc_ol, spp = "Fd", edatope = "C4")



######################## Expanded versions for all edatopes/species/time periods ##########################################################

# Set up
edatopes <- c("C4", "B2", "D6")
base_dir <- "C:/Vanessa/GitHub/SpeciesSummariesFork/Draft_Species_Summaries/Figures"

# 2-panel map set up
bc_ol <- vect("data-raw/data_tables/bc_outline.gpkg")
periods  <- c("2001_2020", "2021_2040", "2041_2060", "2061_2080", "2081_2100")

#Set the species 
# "Pl","Sx","Fd","Cw","Hw","Py","Bl","At","Ac","Ep","Yc","Pw","Ss","Bg","Lw","Mb","Ba"

spp <- "Py"   # Change this for species

## species bubble plots
for (e in edatopes) {
  png(
    filename = file.path(base_dir, spp, paste0(spp, "_", e, "_bub.png")),
    width = 1800, height = 1200, res = 300
  )
  spp_bubbleplot(perexp,
                 period = "2041_2060",
                 scenario = "ssp245",
                 species.focal = spp,
                 edatope = e)
  dev.off()
}

## alluvial plots
for (e in edatopes) {
  png(
    filename = file.path(base_dir, spp, paste0(spp, "_", e, "_alluv.png")),
    width = 2100, height = 1500, res = 300
  )
  g <- plot_spparea(con, spp = spp, edatope = e, fractional = FALSE, save_png = FALSE)
  print(g)
  dev.off()
}


## 2 panel suitability change maps
for (e in edatopes) {
  for (p in periods) {
    
    png(
      filename = file.path(
        base_dir, spp,
        paste0("Two_Panel.", spp, ".", e, ".", p, ".png")
      ),
      width = 1950, height = 1500, res = 300
    )
    
    plot_SuitabilityChangeMap(
      con,
      bgc_template = bgc_template,
      outline = bc_ol,
      spp = spp,
      edatope = e,
      period = p,
      three_panel = FALSE,
      save_png = FALSE
    )
    
    dev.off()
  }
}

##reference suit maps
## USE MIA'S R CODE TO MAKE REF SUIT MAPS

############################################
