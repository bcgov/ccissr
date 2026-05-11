# CCISS_TR_Figures_Bubbleplot_BGC
# 6-panel maps and bubbleplot of BGC change
# Colin Mahony colin.mahony@gov.bc.ca

library(terra)
library(duckdb)
library(data.table)
library(climr)

data("zones_colours_ref")

con <- dbCon_cciss("./bc_2km.duckdb", threads = 8)

setwd("C:/Users/CMAHONY/GitHub/CCISS_repos/ccissr") 

dir <- "" # local copy
dem2 <- rast("dem_BC2kmGrid.tif")
bgcs <- vect(paste0(dir, "BEC13Draft_Simplified.gpkg"))
bgc_template <- make_bgc_template(dem2, bgcs)
dbPopulate(con, bgc_template)

bc_ol <- vect("data-raw/data_tables/bc_outline.gpkg")


zones.bc <- c("BG", "BWBS", "CDF", "CWH", "ESSF", "ICH", "IDF", "MH", "MS", "PP", "SBPS", "SBS", "SWB")

ssp <- "ssp245"
period <- "2041_2060"

ssp.names <- c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")
period.names <- c("2010", "2030", "2050", "2070", "2090")

dir <- "C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/CCISS_Manual/Figures/"
png(filename=paste0(dir, "CCISS.TR.BGCprojections.", period,".png"), type="cairo", units="in", width=6.5, height=6.5, pointsize=10, res=300)

x1 <- c(0.00, 0.30, 0.60, 0.65, 0.65, 0.10)
x2 <- c(0.40, 0.70, 1.00, 1.00, 1.00, 0.63)
y1 <- c(0.65, 0.65, 0.65, 0.325, 0.00, 0.075)
y2 <- c(0.96, 0.96, 0.96, 0.60, 0.275, 0.625)

X <- copy(bgc_template$bgc_rast)

#=============================
## Base plot
par(mar=c(0,0,0,0), bg="white")
plot(0, col="white", xaxt="n", yaxt="n", xlab="", ylab="")
box(col="white")

#=============================
## Mapped BGC zones

par(mar=c(0.1,0.1,1.5,0.1), mgp=c(2,0.25,0))

# par(fig=c(x1[1], x2[1], y1[1], y2[1]), mar=c(0.1,0.1,3,0.1), mgp=c(2,0.25,0), new=TRUE)
par(plt = as.vector(rbind(x1, x2, y1, y2)[,1]), new=TRUE)

dat <- dbGetQuery(con, "select * from bgc_points")

bgc_map(X, dat, 
        boundary = bc_ol,
        legend = TRUE, 
        title = paste("(", letters[1], ") BGC zone map (1961-1990)", sep=""),
        add=TRUE 
)

# mtext(paste("(", letters[1], ") BGC zone map (1961-1990)", sep=""), side=3, line=0.2, adj=0.05, cex=0.8, font=2)

# box()


#=============================
# Recent Period
par(plt = as.vector(rbind(x1, x2, y1, y2)[,2]), new=TRUE)

dat <- dbGetQuery(con, "select * from bgc_raw where period = '2001_2020_obs'") # TODO dummy query until we create obs period bgc preds. 

bgc_map(X, dat[, c("cellnum", "bgc_pred")], 
        boundary = bc_ol,
        add=TRUE,
        title = paste("(", letters[2], ") Observed climate (2001-2020)", sep="")
)

# mtext(paste("(", letters[2], ") Observed climate (2001-2020)", sep=""), side=3, line=0.2, adj=0.05, cex=0.8, font=2)

# exotic.pct <- round(as.numeric(formatC(signif(exotic/length(zone.pred)*100,digits=3), digits=3,format="fg", flag="#")),2)
# legend("topright", cex=0.8, title="Exotic zones", legend=paste(names(exotic), " (", exotic.pct, "%)", sep=""), fill=alpha(ColScheme[as.numeric(factor(names(exotic), zone))], 1), bty="n")


#=============================
## Individual GCM runs

gcms <- list_gcms()[c(5,1,7)]
runs <- c("r4i1p1f1", "r1i1p1f1", "r2i1p3f1")

for(i in 1:3){
  dat <- dbGetQuery(
    con,
    "
    SELECT *
    FROM bgc_raw
    WHERE ssp    = ?
      AND gcm    = ?
      AND run    = ?
      AND period = ?
    ",
    params = list(ssp, gcms[i], runs[i], period)
  )
  
  par(plt = as.vector(rbind(x1, x2, y1, y2)[,i+2]), new = TRUE)
  bgc_map(X, dat[, c("cellnum", "bgc_pred")], 
          boundary = bc_ol,
          label_exotic = 400,
          q_exotic = 0.7,
          add=TRUE,
          title = paste("(", letters[3:5][i], ") ", gcms[i], " (", period.names[which(list_gcm_periods()==period)], ", " , ssp.names[which(list_ssps()==ssp)], ")", sep="")
          )
  # mtext(paste("(", letters[3:5][i], ") ", gcms[i], " (", period.names[which(list_gcm_periods()==period)], ", " , ssp.names[which(list_ssps()==ssp)], ")", sep=""), side=3, line=0.2, adj=0.05, cex=0.8, font=2)
}


# ------------------------
# bubble plot
bgc_perexp <- bgc_persist_expand(con, by_zone = TRUE)

bgc_bubbleplot(bgc_perexp, period = period, plt = as.vector(rbind(x1, x2, y1, y2)[,6]), mar = c(6,10,0.1,0.1))

mtext("(f)", side=3, line=-1.5, adj = 0.025)

dev.off()




# ##=================================
# ###
# 
# # examples of how to query the database
# dbListTables(con)
# dbGetQuery(con, "select*from bgc_raw limit 10")
# dbGetQuery(con, "select*from bgc_points limit 10")
# 
# dat <- dbGetQuery(con, "select * from bgc_raw where ssp = 'ssp245' and gcm = 'ACCESS-ESM1-5' and run = 'ensembleMean' and period = '2041_2060'")
# dat <- dbGetQuery(con, "select * from bgc_raw where ssp = 'ssp245' and gcm = 'GISS-E2-1-G' and period = '2041_2060'")
