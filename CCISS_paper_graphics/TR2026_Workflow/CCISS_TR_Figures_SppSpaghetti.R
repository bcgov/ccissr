# CCISS_TR_Figures_SppSpaghetti
# 3-panel plot of mean suitabilty change for all species
# Colin Mahony colin.mahony@gov.bc.ca

library(terra)
library(duckdb)
library(data.table)
library(climr)

edatopes<- c("B2", "C4", "D6")
edatope.name <- c("Subxeric-poor", "Mesic-medium", "Hygric-rich")
ssp.names <- c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")
period.names <- c("2010", "2030", "2050", "2070", "2090")

con <- dbCon_cciss("./bc_2km.duckdb", threads = 8)

setwd("C:/Users/CMAHONY/GitHub/CCISS_repos/ccissr") 

dir <- "" # local copy
dem2 <- rast("dem_BC2kmGrid.tif")
bgcs <- vect(paste0(dir, "BEC13Draft_Simplified.gpkg"))
bgc_template <- make_bgc_template(dem2, bgcs)
dbPopulate(con, bgc_template)

# compile the species list
data("S1", package = "ccissr")
data("T1", package = "ccissr")
spps <- unique(S1$spp)
spps.candidate <- T1$TreeCode[-which(T1$Exclude=="x" | T1$TreeCode == "X")]
spps <- spps[which(spps%in%spps.candidate)]
spps.native <- spps[which(spps%in%T1[Native=="N", TreeCode])]

# summarized suitaiblity projections for all species
sa <- spp_suit_area(con, spp_list = spps)
setkey(sa, ssp, gcm, run, period)

# unique list of runperiods
runperiods <- unique(sa[, .(ssp, gcm, run, period, MAT_diff)])
setkey(runperiods, ssp, gcm, run, period)

spp.type.names <- c("Boreal", "Temperate", "Mesothermal", "Subalpine", "Exotic")

# Calculate total number of grid cells for calculation of species occupation as a percentage
bgc_points <- values(bgc_template$bgc_rast)
bgc_points <- bgc_template$ids$bgc[match(bgc_points, bgc_template$ids$bgc_id)]
bgc_points <- bgc_points[!is.na(bgc_points)]
BGCs_notin_THLB <- fread("C:/Users/CMAHONY/GitHub/CCISS_repos/ccissr/data-raw/data_tables/BGCs_notin_THLB.csv")
BGCs_notin_THLB <- BGCs_notin_THLB$BGC[which(BGCs_notin_THLB$Exlude=="x")]
totalarea <- length(bgc_points[-which(bgc_points %in% BGCs_notin_THLB)])*(res(dem2)[1]/1000)^2



dir <- "C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/CCISS_Manual/Figures/"
png(filename=paste0(dir, "CCISS.TR.SppSpaghetti.png"), type="cairo", units="in", width=6.5, height=5, pointsize=10, res=300)

mat <- matrix(c(1,2,3,4, 6, 5,5,5),2, byrow=T)   #define the plotting order
layout(mat, widths=c(0.225,1,1,1), heights=c(1, 0.04))   #set up the multipanel plot

par(mar=c(0,0,0,0))
plot(1, type="n", axes=F, xlab="", ylab="")  
text(0.75,1,"Suitable range (% of edatope area in BC)", srt=90, cex=1.2)

x <- c(0,runperiods[, MAT_diff])

for(edatope_sel in edatopes){
  for(spp_sel in spps){
    y <- c(sa[spp == spp_sel & Edatopic== edatope_sel, MappedSuit][1], sa[spp == spp_sel & Edatopic == edatope_sel,.(ssp, gcm, run, period, Proj_Area)][runperiods, on = .(ssp, gcm, run, period), Proj_Area])
    y[is.na(y)] <- 0
    y <- y/totalarea 
    l <- loess(y[order(x)]~x[order(x)])
    assign(paste("line",spp_sel, sep="."), predict(l, seq(0,max(x, na.rm = T), 0.01)))
  }
  
  transform=T
  # for(transform in c(T, F)){
  
  par(mar=c(3.85,0,0,0.2), mgp=c(4, 0.2, 0))
  ylim=if(transform==T) c(-2.05,0.01) else c(0,1)
  plot(0, xlim=c(-3.25,6), ylim=ylim, yaxs="i", xaxs="i", col="white", xaxt="n", yaxt="n", 
       xlab="", 
       ylab="")
  # if(edatope_sel==edatopes[2]){
  #   par(xpd=T)
  #   title(xlab=list(bquote(BC~mean~temperature~change~relative~to~"1961-90"~"("*degree*C*")"), cex=1.2))
  #   par(xpd=F)
  # }
  axis(1, at=0:5, labels = 0:5, tck=0)
  if(edatope_sel==edatopes[1]){
    par(mgp=c(1, 0.2, 0))
    par(xpd=T)
    y.labels <- c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5)
    axis(2, lty=0, at=if(transform==T) log10(y.labels) else seq(0,1,0.1), labels = if(transform==T) paste(y.labels*100, "%", sep="") else format(seq(0,1,0.1), scientific = FALSE, big.mark=","), tck=0, las=2)
    par(xpd=F)
  }
  # rect(-9,0,0, 60000, col="lightgray", border=F)
  # rect(max(x),0,9, ylim[2]*1.1, col="lightgray", border=F)
  
  suit.exotic.final <- vector()
  for(spp_sel in spps[-which(spps%in%spps.native)]){
    line <- get(paste("line",spp_sel, sep="."))
    if(transform==T) line[line<1] <- 1
    if(transform==T) line <- log10(line)
    suit.exotic.final[which(spps[-which(spps%in%spps.native)]==spp_sel)] <- line[length(line)]
  }
  
  suit.native.initial <- vector()
  for(spp_sel in spps[which(spps%in%spps.native)]){
    suit.native.initial[which(spps[which(spps%in%spps.native)]==spp_sel)] <- sa[spp == spp_sel & Edatopic== edatope_sel, MappedSuit][1]
  }
  
  spplist <- spps[which(spps%in%spps.native)][order(suit.native.initial)]
  
  # #Color scheme for species groups
  data("zones_colours_ref")
  boreal <- c("Pl", "Sx", "Sb", "At", "Ep", "Pj", "Acb", "Sw")
  temperate <- c("Fd", "Lw", "Pw", "Py", "Bg", "Act", "Ac")
  mesothermal <- c("Hw", "Cw", "Ba", "Ss", "Dr", "Mb")
  subalpine <- c("Hm", "Yc", "Bl", "Ba")
  ColScheme <- rep(NA, length(spplist))
  ColScheme[which(spplist%in%boreal)] <- as.character(zones_colours_ref$colour[which(zones_colours_ref$classification=="SBS")])
  ColScheme[which(spplist%in%temperate)] <- as.character(zones_colours_ref$colour[which(zones_colours_ref$classification=="IDF")])
  ColScheme[which(spplist%in%mesothermal)] <- as.character(zones_colours_ref$colour[which(zones_colours_ref$classification=="CWH")])
  ColScheme[which(spplist%in%subalpine)] <- as.character(zones_colours_ref$colour[which(zones_colours_ref$classification=="MS")])
  
  if(edatope_sel==edatopes[2]){
    text(-3.2, ylim[1]+0.02, "Boreal species", cex=1.1, srt=90, font=2, pos=4, col=unique(ColScheme[which(spplist%in%boreal)]))
    text(-2.4, ylim[1]+0.02, "Temperate species", cex=1.1, srt=90, font=2, pos=4, col=unique(ColScheme[which(spplist%in%temperate)]))
    text(-1.6, ylim[1]+0.02, "Mesothermal species", cex=1.1, srt=90, font=2, pos=4, col=unique(ColScheme[which(spplist%in%mesothermal)]))
    text(-0.8, ylim[1]+0.02, "Subalpine species", cex=1.1, srt=90, font=2, pos=4, col=unique(ColScheme[which(spplist%in%subalpine)]))
  }
  
  for(spp_sel in spplist){
    i <- which(spplist==spp_sel)
    line <- get(paste("line",spp_sel, sep="."))
    if(transform==T) line[line<10^(-5)] <- 10^(-5)
    if(transform==T) line <- log10(line)
    if(line[1]> if(transform==T) ylim[1] else 100){
      lines(seq(0,max(x), 0.01), line, col=ColScheme[i], lwd=2)
      position <- rep(0:3, times=100)
      text(0-position[i]*0.8, line[1], spp_sel, pos=2, col=ColScheme[i], font=2, cex=0.9, offset=0.1)
      lines(c(0-position[i]*0.8,0), rep(line[1],2), col=ColScheme[i], lty=2)
    }
  }
  
  spplist <- spps[-which(spps%in%spps.native)][order(suit.exotic.final)]
  for(spp_sel in spplist){
    i <- which(spplist==spp_sel)
    line <- get(paste("line",spp_sel, sep="."))
    if(transform==T) line[line<10^(-5)] <- 10^(-5)
    if(transform==T) line <- log10(line)
    if(max(line)> if(transform==T) ylim[1] else 100){
      lines(seq(0,max(x), 0.01), line)
      position <- rep(0, times=100)
      if(which.max(line)>(length(line)-100)){
        text(max(x)+position[i]*0.55, line[length(line)], spp_sel, pos=4, cex=0.9, offset=0.1, font=2)
        lines(c(max(x), max(x)+position[i]*0.55), rep(line[length(line)],2), lty=2, lwd=0.6)
      } else {
        text(seq(0,max(x), 0.01)[which(line==max(line))]+position[i]*0.2, max(line), spp_sel, pos=3, cex=0.9, offset=0.1, font=2)
      }
    }
  }
  rect(0,-10,max(x),ylim[2]*1.1, col=NA, border=T)
  box()
  
  
  
  # boxplot for focal period
  if(edatope_sel==edatopes[2]){
    par(xpd=T)
    ssp.focal <- "ssp245"
    for(i in 1:4){
      ssp.focal <- list_ssps()[c(2,2,2,3)][i]
      period.focal <- list_gcm_periods()[c(1,3,5,5)][i]
      x.focal <- runperiods[ssp==ssp.focal & period==period.focal, MAT_diff]
      position <- ylim[1] - diff(ylim)/40 - diff(ylim)/50*i
      boxplot(x.focal, add=T, col=c("dodgerblue", "red")[which(list_ssps()[2:3]==ssp.focal)], horizontal=TRUE, axes=FALSE, range=0, at=position, boxwex = diff(ylim)/50)
      text(if(ssp.focal=="ssp370" & period.focal==list_gcm_periods()[5]) min(x.focal) else max(x.focal), 
           position, 
           paste(ssp.names[which(list_ssps()==ssp.focal)], ", ", period.names[which(list_gcm_periods()==period.focal)], sep=""), 
           pos=if(ssp.focal=="ssp370" & period.focal==list_gcm_periods()[5]) 2 else 4, cex=0.9)
    }
    par(xpd=F)
  }
  
  
  mtext(paste("(", letters[which(edatopes==edatope_sel)],") ", edatope_sel, " edatope", sep=""), side=3, line=-1.5, adj=0.58, cex=0.8, font=2)
  mtext("Native", side=3, line=-1.5, adj=0.025, cex=0.8, font=2)
  mtext("Non-native", side=4, line=-1.5, adj=0.975, cex=0.8, font=2)
  
  # }
  print(edatope_sel)

}
par(mar=c(0,0,0,0))
plot(1, type="n", axes=F, xlab="", ylab="")  
text(1,1,bquote(BC~mean~temperature~change~relative~to~"1961-90"~"("*degree*C*")"), srt=0, cex=1.2)

dev.off()



