# 1-panel bubbleplot of subzone/variant
# Colin Mahony colin.mahony@gov.bc.ca

library(terra)
library(duckdb)
library(data.table)
library(climr)

con <- dbCon_cciss("./bc_2km.duckdb", threads = 8)

setwd("C:/Users/CMAHONY/GitHub/CCISS_repos/ccissr") 

dir <- "" # local copy
dem2 <- rast("dem_BC2kmGrid.tif")
bgcs <- vect(paste0(dir, "BEC13Draft_Simplified.gpkg"))
bgc_template <- make_bgc_template(dem2, bgcs)
dbPopulate(con, bgc_template)

ssp <- "ssp245"
period.focal <- "2021_2040"

ssp.names <- c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")
period.names <- c("2001-2020", "2021-2040", "2041-2060", "2061-2080", "2081-2100")

  
# ------------------------
# bubble plot

scenario = "ssp245"
xlab = "Climate analog persistence within historical range"
ylab = "Climate analog expansion beyond historical range"
unit.focal = NULL
focal.color = "lightskyblue2"
xlim = c(-8, 0.3)
ylim = c(-6,3)
xlabels = TRUE
ylabels = TRUE
mar = c(3,4,0.1,0.1) 
plt = NULL

bgc_perexp <- bgc_persist_expand(con, by_zone = FALSE)
persist_expand <- na.omit(bgc_perexp, col = c("Persistance","Expansion"))

units <- unique(persist_expand$bgc)

persist_expand[Expansion<2^(ylim[1]-1), Expansion := 2^(ylim[1])]
persist_expand[, Expansion := log2(Expansion)]
persist_expand[Persistance<2^(xlim[1]-1), Persistance := 2^(xlim[1])]
persist_expand[, Persistance := log2(Persistance)]


dir <- "C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/CCISS_Manual/Figures/"
png(filename=paste0(dir, "CCISS.TR.SZbubbleplot.png"), type="cairo", units="in", width=6.5, height=8, pointsize=10, res=300)

layout(matrix(c(
  1, 2, 3, 4,
  5, 5, 5, 5
), nrow = 2, byrow = TRUE),
heights = c(1, 4))

# Outer margins for shared labels
par(oma = c(2.5, 3.5, 1, 0), mar = c(.21, .21, .21, .21), mgp=c(1.5, 0.25, 0), cex=1)

# plots
it.focal <- which(list_gcm_periods()==period.focal)
its <- c(1:4, it.focal)
for(i in 1:5) {

  period_sel <- list_gcm_periods()[its][i]
    
  plot(0, xlim=xlim, ylim=ylim, col="white", xaxt="n", yaxt="n", xlab= if(i == 5) xlab else "", ylab="", xaxs="i")
  
  lines(log2(1-iso), log2(iso), lty=2, lwd=2, col="darkgray")
  
  if(i == 1) axis(2,at=seq(ylim[1], ylim[2], 2), labels=paste(round(2^(seq(ylim[1], ylim[2],2))*100),"%", sep=""), las=2, tck=0)
  
  for(unit in units){
    col.focal <- WNA_BGCs$ZoneColour[which(WNA_BGCs$BGC==unit)] 
    col.focal2 <- "black" 
    x <- persist_expand[(ssp == scenario | is.na(ssp)) & period == period_sel & bgc_pred == unit, Persistance]
    y <- persist_expand[(ssp == scenario | is.na(ssp)) & period == period_sel & bgc_pred == unit, Expansion]
    
    # points(x,y)
    xuse <- mean(x) 
    yuse <- mean(y)
    xuse[xuse > (xlim[2] - 0.01)] <- xlim[2] - 0.01
    yuse[yuse > (ylim[2] - 0.01)] <- ylim[2] - 0.01
    points(xuse,yuse, pch=21, bg=col.focal, cex= 1, col=col.focal2)
    if(i==5) text(xuse,yuse, unit, cex= 0.6, font=1, col=col.focal2, pos=4, offset=0.3)
  }
  
  mtext(paste0("(", letters[i], ") ", period.names[list_gcm_periods()==period_sel]), side=3, line= if(i==5) -1.5 else 0.25, adj = 0.025)
  
  box()
}

axis(1,at=seq(xlim[1], xlim[2]), labels=paste(round(2^(seq(xlim[1], xlim[2]))*100),"%", sep=""), tck=0)
axis(2,at=seq(ylim[1], ylim[2]), labels=paste(round(2^(seq(ylim[1], ylim[2]))*100),"%", sep=""), las=2, tck=0)
mtext(xlab, side = 1, outer = TRUE, line = 1.25)
mtext(ylab, side = 2, outer = TRUE, line = 2.5)



dev.off()
