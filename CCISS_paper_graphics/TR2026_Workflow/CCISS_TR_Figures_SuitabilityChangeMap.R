# CCISS_TR_Figures_SuitabilityChangeMap
# 2-panel map of suitabilty change for the CCISS technical report
# Colin Mahony colin.mahony@gov.bc.ca

library(ccissr)
library(data.table)
library(terra)
library(duckdb)

con <- dbCon_cciss("./bc_2km.duckdb", threads = 8)

setwd("C:/Users/CMAHONY/GitHub/CCISS_repos/ccissr") 
bc_ol <- vect("data-raw/data_tables/bc_outline.gpkg")

dem2 <- rast("dem_BC2kmGrid.tif")

dir <- "" # local copy
bgcs <- vect(paste0(dir, "BEC13Draft_Simplified.gpkg"))
bgc_template <- make_bgc_template(dem2, bgcs)

# compile the species list
data("S1", package = "ccissr")
data("T1", package = "ccissr")
spps <- unique(S1$spp)
spps.candidate <- T1$TreeCode[-which(T1$Exclude=="x" | T1$TreeCode == "X")]
spps <- spps[which(spps%in%spps.candidate)]
spps.native <- spps[which(spps%in%T1[Native=="N", TreeCode])]

# populate any missing species/periods
summarise_preds(con)
siteseries_preds(con)
cciss_suitability(con, spps)

edatopes<- c("B2", "C4", "D6")

edatope_sel = "C4"
for(edatope_sel in edatopes){
  
  spp_sel = "Fd"
  for(spp_sel in spps){
    
    period_sel <- "2001_2020"
    period_sel <- "2021_2040"
    for(period_sel in list_gcm_periods()){
      
      dir <- "C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/CCISS_Manual/Figures/SuitChangeMaps/"
      png(filename=paste0(dir, paste("CCISS.TR.SuitabilityChangeMap", edatope_sel, spp_sel, period_sel,"png", sep=".")), type="cairo", units="in", width=6.5, height=5, pointsize=12, res=300)
      
      plot_SuitabilityChangeMap(con, spp = spp_sel, edatope = edatope_sel, period = period_sel, bgc_template = bgc_template, outline = bc_ol, three_panel = FALSE, save_png = FALSE)
      
      dev.off()
      
      print(period_sel)
    }
    print(spp_sel)
  }
  print(edatope_sel)
}

# # Troubleshooting scripts
# dbGetQuery(con, "select * from cciss_res where Spp = 'Fd' and FuturePeriod = '2001_2020' limit 10")
# dbExecute(con, "drop table cciss_res") # deletes the cciss_res table so that it can be built from scratch. 
