
# checking whether MAT_diff from ccissr::spp_suit_area is correct [it isn't]
# run in the ccissr package environment


dem2 <- rast("dem_BC2kmGrid.tif")
dem <- aggregate(dem2, 5)
dem_table <- climr::dem_to_table(dem)

# -------------------------------
# climr data

data_climr <- downscale(dem_table, 
                        gcms = list_gcms()[c(1,4,5,6,7,10,11,12)], 
                        gcm_periods = list_gcm_periods(), 
                        ssps = list_ssps()[1:3], 
                        max_run = 4L, 
                        vars = "MAT"
)

data_climr[, MAT_diff := MAT - MAT[1L], by = id]

mean_climr <- data_climr[ , .(
  MAT = mean(MAT, na.rm = TRUE),
  MAT_diff = mean(MAT_diff, na.rm = TRUE)
),
by = .(GCM, SSP, PERIOD, RUN)
][-1,]

setkey(mean_climr, GCM, SSP, RUN, PERIOD)
mean_climr

hist(mean_climr$MAT_diff)

# -------------------------------
# ccissr data

con <- dbCon_cciss("./bc_2km.duckdb", threads = 8)

sa <- spp_suit_area(con, spp_list = "Fd")
mean_ccissr <- unique(sa[, .(ssp, gcm, run, period, MAT_diff)])
setnames(mean_ccissr,
         old = c("gcm", "ssp", "run", "period"),
         new = c("GCM", "SSP", "RUN", "PERIOD"))
setkey(mean_ccissr, GCM, SSP, RUN, PERIOD)


# -------------------------------
# difference

diff <- mean_climr[
  mean_ccissr,
  MAT_diff - i.MAT_diff
]
hist(diff)
