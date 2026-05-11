# CCISS_TR_Figures_Bubbleplot
# Three-panel bubbleplot figure for the CCISS technical report
# Colin Mahony colin.mahony@gov.bc.ca

library(ccissr)
library(data.table)
library(terra)

edatopes<- c("B2", "C4", "D6")
edatope.name <- c("Subxeric-poor", "Mesic-medium", "Hygric-rich")


dem2 <- rast("dem_BC2kmGrid.tif")

dir <- "" # local copy
bgcs <- vect(paste0(dir, "BEC13Draft_Simplified.gpkg"))
bgc_template <- make_bgc_template(dem2, bgcs)

con <- dbCon_cciss("./bc_2km.duckdb", threads = 8)
dbPopulate(con, bgc_template)

perexp <- spp_persist_expand(con, spp_list = c("Pl", "Fd", "Cw", "Sx", "At", "Py", "Ba", "Bl", "Bg", "Yc", "Hm", "Lw", "Hw","Dr", "Ep", "Act", "Sb", "Mb", "Ss"), fractional = TRUE)

period_sel <- "2021_2040"
for(period_sel in list_gcm_periods()){
  
  dir <- "C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/CCISS_Manual/Figures/"
  png(filename=paste0(dir, "CCISS.TR.BubblePlot.THLB.", period_sel,".png"), type="cairo", units="in", width=6.5, height=4.25, pointsize=9, res=400)
  mat <- matrix(c(1,1,2,1,1,3,5,5,5,4,4,4),4, byrow=T)   #define the plotting order
  layout(mat, widths=c(1,1,1), heights=c(1,1,0.1,0.1))   #set up the multipanel plot
  
    for(edatope_sel in edatopes[c(2,1,3)]){
    
    args <- list(
      perexp,
      period = period_sel,
      edatope = edatope_sel,
      mar = c(0.1, 4, 0.1, 0.1),
      xlab = "" 
    )
    
    if (edatope_sel != "C4") {
      args$mar <- c(0.1, 0.1, 0.1, 0.1)
      args$ylab <- ""
      args$xlabels <- FALSE
      args$ylabels <- FALSE
      annotate = FALSE
    }
    
    if (edatope_sel == "D6") {
      args$xlabels <- TRUE
    }
    
    do.call(spp_bubbleplot, args) 
    
    panel <- paste("(", letters[which(edatopes[c(2,1,3)]==edatope_sel)],")", sep="")
    mtext(paste(panel," ", edatope_sel," edatope", " (", edatope.name[which(edatopes==edatope_sel)], " sites)", sep=""), side=3, line=-1.25, adj= if(edatope_sel=="C4") 0.025 else 0.075, cex=1, font=1)
    
  }
  
  par(mar=c(0,0,0,0))
  plot(1, type="n", axes=F, xlab="", ylab="")  
  text(1,1, "Persistence within historically suitable range",cex=1.2)  
  
  dev.off()
  print(period_sel)
}

