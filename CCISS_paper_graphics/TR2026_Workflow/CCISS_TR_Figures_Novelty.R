# CCISS_TR_Figures_Bubbleplot_BGC
# 6-panel maps and bubbleplot of BGC change
# Colin Mahony colin.mahony@gov.bc.ca

library(terra)
library(duckdb)
library(data.table)
library(climr)
library(ccissr)

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
period.names <- c("2001-2020", "2021-2040", "2041-2060", "2061-2080", "2081-2100")

# Color Scheme for sigma novelty
breakseq <- c(0,4,8)
breakpoints <- c(seq(breakseq[1], breakseq[3], 0.01),199); length(breakpoints)
ColScheme <- c(colorRampPalette(c("gray90", "gray50", "#FFF200", "#CD0000", "black"))(length(breakpoints)))

dir <- "C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/CCISS_Manual/Figures/"
png(filename=paste0(dir, "CCISS.TR.Novelty.", period,".png"), type="cairo", units="in", width=6.5, height=5.4, pointsize=10, res=300)

# x1 <- c(0.00, 0.30, 0.60, 0.65, 0.65, 0.10)
# x2 <- c(0.40, 0.70, 1.00, 1.00, 1.00, 0.63)
# y1 <- c(0.65, 0.65, 0.65, 0.325, 0.00, 0.075)
# y2 <- c(0.96, 0.96, 0.96, 0.60, 0.275, 0.625)

x1 <- c(0.00, 0.425, 0.00, 0.2375, 0.4675, 0.70)
x2 <- c(0.575, 1.00, 0.30, 0.5375, 0.7675, 1.00)
y1 <- c(0.36, 0.36, 0.00, 0.00, 0.00, 0.00)
y2 <- c(0.96, 0.96, 0.30, 0.30, 0.30, 0.30)

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

dat <- dbGetQuery(con, "select * from bgc_raw where ssp = 'ssp245' and gcm = 'EC-Earth3' and run = 'ensembleMean' and period = '2041_2060'")

bgc_map(X, dat[, c("cellnum", "bgc_pred")], 
        # boundary = bc_ol,
        label_exotic = 400,
        q_exotic = 0.7,
        legend = TRUE, 
        title = paste("(", letters[1], ") Ensemble BGC Projection, ", period.names[which(list_gcm_periods()==period)], sep=""),
        add=TRUE 
)


# mtext(paste("(", letters[1], ") BGC zone map (1961-1990)", sep=""), side=3, line=0.2, adj=0.05, cex=0.8, font=2)

# box()


#=============================
# Novelty for selected Period
par(plt = as.vector(rbind(x1, x2, y1, y2)[,2]), new=TRUE)
title <- paste("(", letters[2], ") Ensemble climatic novelty, ", period.names[which(list_gcm_periods()==period)], sep="")
novelty <- rast(paste0("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/CCISS_Manual/Data/ClimaticNovelty_Provincial_SZ_Ensemble_", period, ".tif"))  
novelty <- project(novelty, X)
image(novelty, axes=F, col=ColScheme, main = title , adj = 0.05, cex.main = 0.85, font.main = 1)
# plot(novelty, legend=F, axes=F, col=ColScheme, add=T)
# plot(bc_ol, add=T, lwd=0.4)

legend_ramp(
  r = X,
  title = "Sigma Novelty",
  ColScheme = ColScheme,
  breaks = breakpoints[-length(breakpoints)],
  pos = c(0.8, 0.85, 0.45, 0.95),
  log = NULL,
  log.relative = FALSE,
  horizontal = FALSE,
  title.height = 4,
  margin = 0.02
)

#=============================
## Novelty for other periods

index <- 3 #plot counter
for(i in c(1:5)[-which(list_gcm_periods() == period)]){
  novelty <- rast(paste0("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/CCISS_Manual/Data/ClimaticNovelty_Provincial_SZ_Ensemble_", list_gcm_periods()[i], ".tif"))  
  novelty <- project(novelty, X)

  par(plt = as.vector(rbind(x1, x2, y1, y2)[,index]), new = TRUE)
  title <- paste("(", letters[index], ") ", period.names[i], sep="")
  image(novelty, axes=F, col=ColScheme, main = title , adj = 0.05, cex.main = 0.85, font.main = 1)
  # plot(bc_ol, add=T, lwd=0.4)

  index <- index + 1
  
}


dev.off()




# ##=================================
# ###
# 
# # examples of how to query the database
# dbListTables(con)
# dbGetQuery(con, "select*from bgc_raw limit 1000")
# dbGetQuery(con, "select*from bgc_points limit 10")
# dbGetQuery(con, "select*from bgc_raw where ssp = 'ssp245' and gcm = 'ACCESS-ESM1-5' and run = 'ensembleMean' and period = '2041_2060'")
# dbGetQuery(con, "select*from bgc_raw where period = '2041_2060'")
# dbGetQuery(con, "
#   SELECT DISTINCT run
#   FROM bgc_raw
# ")
# 
# dat <- dbGetQuery(con, "select * from bgc_raw where ssp = 'ssp245' and gcm = 'ACCESS-ESM1-5' and run = 'ensembleMean' and period = '2041_2060'")
# dat <- dbGetQuery(con, "select * from bgc_raw where ssp = 'ssp245' and gcm = 'GISS-E2-1-G' and period = '2041_2060'")
