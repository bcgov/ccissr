### Kiri Daust, 2025

#' Create a spaghetti plot for a single species over multiple edatopes
#' @param suit_area data.table. Usually created with `spp_suit_area`
#' @param species Character. Single species to create plot for
#' @param use_MAT Use change in mean annual temperature as x-axis instead of time? Default to `FALSE`
#' @return NULL. Creates plot
#' @import data.table
#' @import ggplot2
#' @importFrom stinepack stinterp
#' @export
spp_spaghettiplot <- function(suit_area, species, use_MAT = FALSE) {
  suit_area <- suit_area[spp == species,]
  if(use_MAT) {
    suit_area[,xvar := MAT_diff]
  } else {
    suit_area[,xvar := as.integer(substr(period,1,4))]
  }

  dat_spline <- suit_area[
    , {
      if(use_MAT){
        dt <- rbind(.SD, data.table(xvar = 0, Suit_Prop = 1), fill = TRUE)
      } else {
        dt <- rbind(.SD, data.table(xvar = 2000, Suit_Prop = 1), fill = TRUE)
      }
      o <- order(dt$xvar)
      x <- dt$xvar[o]
      y <- dt$Suit_Prop[o]
      
      # interpolate
      xout <- seq(min(x), max(x), length.out = 100)
      yout <- stinterp(x, y, xout)$y
      
      .(xvar = xout, Suit_Spline = yout, Spline_Num = 1:100)
    },
    by = .(Edatopic, ssp, gcm, run)
  ]
  dat_spline[, Group := interaction(ssp, gcm, run)]
  if(use_MAT){
    mean_spline <- dat_spline[
      , .(Suit_Spline = mean(Suit_Spline, na.rm = T), xvar = mean(xvar, na.rm = T)),
      by = .(Edatopic, Spline_Num)
    ]
    ggplot(dat_spline, aes(x = xvar, y = Suit_Spline, group = Group)) +
      geom_line(alpha = 0.4) +
      geom_line(data = mean_spline, aes(x = xvar, y = Suit_Spline), inherit.aes = FALSE, linewidth = 1.4) +
      facet_wrap(~Edatopic)+
      theme_minimal() +
      labs(y="Proportion of Historic Suitable Area", x="Change in mean annual temperature")
  } else {
    mean_spline <- dat_spline[
      , .(Suit_Spline = mean(Suit_Spline), xvar = mean(xvar)),
      by = .(Edatopic, ssp, Spline_Num)
    ]
    ggplot(dat_spline, aes(x = xvar, y = Suit_Spline, group = Group, colour = ssp)) +
      geom_line(alpha = 0.2) +
      geom_line(data = mean_spline, aes(group = ssp), size = 1.4) +
      facet_wrap(~Edatopic)+
      theme_minimal() +
      labs(y="Proportion of Historic Suitable Area", x="Time Period", colour="SSPs") +
      scale_colour_manual(values = c("#E6AB02", "#D95F02", "#7F0000"))
  }
}


#' Create C.R. Mahony's bubbleplot of species expansion/persistance
#' @param persist_expand data.table. Usually created with `spp_persist_expand`
#' @param period_sel Character. Period to create plot for.
#' @param scenario Character. ssp to use for plot.
#' @param species Character vector of species or "auto" to display all available species. Default "auto".
#' @param edatope Character of desired edatope to use for multiple species plot, or vector of edatopes for single species.
#' @param by  Character. Either "species", or "edatopes". If 
#' @return NULL. Creates plot
#' @import data.table
#' @importFrom car dataEllipse
#' @importFrom RColorBrewer brewer.pal
#' @export
spp_bubbleplot <- function(persist_expand, period_sel, scenario, species = "auto", edatope, by = "species") {
  
  if(species == "auto"){
    spps <- unique(persist_expand$spp)
  } else if (by == "species") {
    spps <- species
  } else {
    spps <- edatope
  }
  #spps <- unique(persist_expand$spp)
  colors = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)][-1]
  colors = colors[-grep("yellow", colors)]
  set.seed(5)
  if(length(spps) > 11){
    sppcolors <- c(brewer.pal(n=12, "Paired")[-11],sample(colors,length(spps) - 11)) # removal of "11" is light yellow, doesn't show up well.
  } else {
    sppcolors <- c(brewer.pal(n=12, "Paired")[-11])
  }
  
  par(mar=c(3,4,0.1,0.1), mgp=c(1.25, 0.25, 0), cex=1.5)
  
  xlim <- c(0, 1.5)
  ylim <- c(-5,3)
  plot(0, xlim=xlim, ylim=ylim, col="white", xaxt="n", yaxt="n", xlab="Persistence within historically feasible range", ylab="")
  axis(1,at=seq(xlim[1], xlim[2], 0.2), labels=paste(seq(xlim[1], xlim[2], 0.2)*100,"%", sep=""), tck=0)
  axis(2,at=seq(ylim[1], ylim[2]), labels=paste(round(2^(seq(ylim[1], ylim[2]))*100),"%", sep=""), las=2, tck=0)
  par(mgp=c(2.75, 0.25, 0))
  title(ylab="Expansion beyond historically feasible range", cex.lab=1)
  iso <- seq(0,1.2, 0.001)
  lines(1-iso, log2(iso), lty=2, lwd=2, col="darkgray")
  #arctext(x = "Growing feasible range", center = c(-1, -28.7), radius = 4.6, start = 0.431*pi , cex = 0.8, stretch = 1.05, col="darkgray", font=2)
  #arctext(x = "Shrinking feasible range", center = c(-1, -29.3), radius = 4.6, start = 0.431*pi , cex = 0.8, stretch = 1.05, col="darkgray", font=2)
  #mtext(paste(edatope.names[which(edatopes==edatope)], " sites", " (", edatope, ")", sep=""), side=3, line=-1.25, adj= if(edatope=="C4") 0.025 else 0.075, cex=0.7, font=1)
  
  
  spp.focal <- "none"
  #spps <- spps[-9]
  for(i in 1:length(spps)){
    if(by == "species") {
      spp_sel <- spps[i]
      eda_sel <- edatope[1]
    } else {
      spp_sel <- species[1]
      eda_sel <- spps[i]
    }
    col.focal <- if(spp.focal=="none") sppcolors[i] else "lightgray"
    col.focal2 <- if(spp.focal=="none") "black" else "darkgray"
    x <- persist_expand[ssp == scenario & period == period_sel & Edatopic == eda_sel & spp == spp_sel, Persistance]
    y <- persist_expand[ssp == scenario & period == period_sel & Edatopic == eda_sel & spp == spp_sel, Expansion]
    y[y<2^(ylim[1])] <- 2^(ylim[1])
    y <- log2(y)
    
    # points(x,y)
    if(length(x)>1 & var(x)>0){
      if(var(y)==0) {
        lines(range(x), range(y), col=col.focal)
      }  else {
        dataEllipse(x, y, levels=0.5, center.pch=21, add=T, col=col.focal, fill=T, lwd=0.5, plot.points=F)
      } 
    }
    points(mean(x),mean(y), pch=21, bg=col.focal, cex=if(spps[i]==spp.focal) 4.5 else 3, col=col.focal2)
    text(mean(x),mean(y), spps[i], cex=if(spps[i]==spp.focal) 1 else 0.7, font=2, col=col.focal2)
  }
  
}


#' Create C.R. Mahony's bubbleplot of bgc expansion/persistance
#' @param persist_expand data.table. Usually created with `bgc_persist_expand`
#' @param period Character. Period to create plot for.
#' @param scenario Character. ssp to use for plot.
#' @return NULL. Creates plot
#' @import data.table
#' @importFrom car dataEllipse
#' @export
bgc_bubbleplot <- function(persist_expand, period, scenario) {
  unit.persistence.focal <- "none"
  persist_expand <- na.omit(persist_expand, col = c("Persistance","Expansion"))
  ColScheme <- rbind(copy(subzones_colours_ref),copy(zones_colours_ref))
  
  units <- unique(persist_expand$bgc)
  period_sel <- period
  
  par(mar=c(3,4,0.1,0.1), mgp=c(1.25, 0.25, 0), cex=1.5)
  
  xlim <- c(0, 1.1)
  ylim <- c(-5,3)
  plot(0, xlim=xlim, ylim=ylim, col="white", xaxt="n", yaxt="n", xlab="Persistence within historical range", ylab="")
  axis(1,at=seq(xlim[1], xlim[2], 0.2), labels=paste(seq(xlim[1], xlim[2], 0.2)*100,"%", sep=""), tck=0)
  axis(2,at=seq(ylim[1], ylim[2]), labels=paste(round(2^(seq(ylim[1], ylim[2]))*100),"%", sep=""), las=2, tck=0)
  par(mgp=c(2.75, 0.25, 0))
  title(ylab="Expansion beyond historical range", cex.lab=1)
  iso <- seq(0,1.2, 0.001)
  lines(1-iso, log2(iso), lty=2, lwd=2, col="darkgray")
  
  for(unit in units){
    col.focal <- if(unit.persistence.focal=="none") ColScheme$colour[which(ColScheme$classification==unit)] else "lightgray"
    col.focal2 <- if(unit.persistence.focal=="none") "black" else "darkgray"
    x <- persist_expand[ssp == scenario & period == period_sel & bgc_pred == unit, Persistance]
    y <- persist_expand[ssp == scenario & period == period_sel & bgc_pred == unit, Expansion]
    y[y<2^(ylim[1])] <- 2^(ylim[1])
    y <- log2(y)
    
    # points(x,y)
    if(length(x)>1 & var(x) > 0){
      if(var(y)==0) {
        lines(range(x, na.rm=T), range(y), col=col.focal)
      } else {
        dataEllipse(x, y, levels=0.5, center.pch=21, add=T, col=col.focal, fill=T, lwd=0.5, plot.points=F)
      } 
    }
    points(mean(x),mean(y), pch=21, bg=col.focal, cex=if(unit==unit.persistence.focal) 4.5 else 3, col=col.focal2)
    text(mean(x),mean(y), unit, cex=if(unit==unit.persistence.focal) 1 else 0.7, font=2, col=col.focal2)
  }
  
}

#' Create alluvial/stacked bar plot of projected area by zone/subzone
#' @param spp Character. Single species to use for plot
#' @param edatope Character. Single edatopic position for plot (e.g., "C4")
#' @param bgc_template List containing SpatRaster of BGCs and id table, or data.table (must have columns `cell`,`BGC`). Usually created using `make_bgc_template`
#' @param fractional Logical. Use fractional (suitability based) values for calculations?
#' @param by_zone Logical. Plot by zone or subzone? Defaults to `TRUE`
#' @param base_folder Base folder to read results from.
#' @return NULL. Creates plot
#' @import data.table duckdb ggplot2
#' @importFrom ggalluvial geom_alluvium
#' @export

plot_spparea <- function(dbCon, spp, edatope, fractional, by_zone = TRUE) {
  cciss_spp <- dbGetQuery(dbCon, sprintf("select * from cciss_res where Spp = '%s' AND Edatope = '%s'", spp, edatope)) |> as.data.table()
  bgc_mapped <- dbGetQuery(dbCon, "select * from bgc_points") |> as.data.table()
  
  if(by_zone) {
    bgc_mapped[, zone := regmatches(bgc, regexpr("^[A-Z]+", bgc))]
  } else {
    bgc_mapped[, zone := bgc]
  }
  
  cciss_spp[bgc_mapped, zone := i.zone, on = c(SiteRef = "cellnum")]
  
  
  if(fractional) {
    cciss_spp[,SppSuit := Newsuit]
    cciss_spp[is.na(SppSuit) | SppSuit > 3.5, SppSuit := 5]
    cciss_spp[, SppSuit := 1 - (SppSuit - 1) / 4]
    
    cciss_spp[,HistSuit := Curr]
    cciss_spp[is.na(HistSuit) | HistSuit > 3.5, HistSuit := 5]
    cciss_spp[, HistSuit := 1 - (HistSuit - 1) / 4]
  } else {
    cciss_spp[,SppSuit := 0]
    cciss_spp[Newsuit <= 3, SppSuit := 1]
    
    cciss_spp[,HistSuit := 0]
    cciss_spp[Curr <= 3, HistSuit := 1]
  }
  
  cciss_sum <- cciss_spp[,.(SppArea = sum(SppSuit)), by = .(zone, FuturePeriod)]
  cciss_hist <- cciss_spp[,.(SppArea = sum(HistSuit)), by = .(zone, FuturePeriod)][FuturePeriod == "2021_2040",][,FuturePeriod := "1961"]
  cciss_sum <- rbind(cciss_sum, cciss_hist)
  cciss_sum[, Year := as.factor(substr(FuturePeriod,1,4))]
  
  yrs <- unique(as.integer(cciss_sum$Year))
  idx <- seq_along(yrs)
  
  midpts <- head(idx, -1) + 0.5
  
  bars <- data.frame(
    xmin = midpts - 0.3,
    xmax = midpts + 0.3,
    ymin = -Inf,
    ymax = Inf
  )
  
  if(by_zone) {
    colScheme <- c(PP = "#ea7200", MH = "#6f2997", SBS = "#2f7bd2", ESSF = "#ae38b8", 
                     CWH = "#488612", BWBS = "#4f54cf", CWF = "#7577e7", IGF = "#77a2eb", 
                     CMX = "#71d29e", BG = "#dd1320", IDF = "#e5d521", MS = "#e44ebc", 
                     SWB = "#a1dbde", CRF = "#af3a13", WJP = "#73330e", ICH = "#1fec26", 
                     CDF = "#edf418", JPW = "#96b3a5", CMA = "#eae1ee", SBPS = "#6edde9", 
                     IMA = "#e3f1fa", GBD = "#4d433f", OW = "#582511", BAFA = "#eee4f1", 
                     MMM = "#FF00FF", MHRF = "#2612dc", MGP = "#f0aeab", FG = "#92696c", 
                     SGP = "#cca261", GO = "#f0a325", SBAP = "#51d5a7", IWF = "#d44273", 
                     BSJP = "#424160", MSSD = "#dac370", MDCH = "#2d0cd4", CVG = "#c9edd3", 
                     SAS = "#92b1b6", CCH = "#7e22ca")
  } else {
    colScheme <- setNames(subzones_colours_ref$colour, subzones_colours_ref$classification)
  }
  
  # cellarea <- (res(bgc_template$bgc_rast)[2]*111)*(res(bgc_template$bgc_rast)[1]*111*cos(mean(ext(bgc_template$bgc_rast)[3:4]) * pi / 180))
  # cciss_sum[, SppArea := SppArea * cellarea]
  
  ggplot(cciss_sum, aes(x = Year, y = SppArea, fill = zone)) +
    geom_alluvium(aes(alluvium = zone), alpha= .9, color = "black") +
    geom_rect(
      data = bars,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "white", alpha = 0.9    # adjust alpha for translucency
    ) +
    theme_classic() +
    scale_fill_manual(values = colScheme) +
    labs(y="Species Suitable Area (Km^2)",x="Time period")
  
}

#' Create C.R. Mahony's two-panel map plot of historic suitability, suitability change, and a boxplot of change by zone
#' @description
#' Note that currently, this function only works correctly if the preceeding analysis has been done using an Albers grid.
#' @param dbCon duckdb database connection
#' @param spp Character. Species code to create plot for
#' @param edatope Character. Edatope code to create plot for (e.g., "C4")
#' @param period Character. Period to create plot for.
#' @param bgc_template List containing SpatRaster of BGCs and id table. Usually created using `make_bgc_template`
#' @param outline SpatVector of aoi boundary.
#' @param save_png Logical. Save plot to png? Default `TRUE`. If `FALSE`, creates plot on default plotting device.
#' @return NULL. Creates plot
#' @import data.table duckdb terra
#' @export
plot_2panel <- function(dbCon, spp, edatope, period, bgc_template, outline, three_panel = FALSE, save_png = TRUE) {
  
  zoneScheme <- c(PP = "#ea7200", MH = "#6f2997", SBS = "#2f7bd2", ESSF = "#ae38b8", 
                 CWH = "#488612", BWBS = "#4f54cf", CWF = "#7577e7", IGF = "#77a2eb", 
                 CMX = "#71d29e", BG = "#dd1320", IDF = "#e5d521", MS = "#e44ebc", 
                 SWB = "#a1dbde", CRF = "#af3a13", WJP = "#73330e", ICH = "#1fec26", 
                 CDF = "#edf418", JPW = "#96b3a5", CMA = "#eae1ee", SBPS = "#6edde9", 
                 IMA = "#e3f1fa", GBD = "#4d433f", OW = "#582511", BAFA = "#eee4f1", 
                 MMM = "#FF00FF", MHRF = "#2612dc", MGP = "#f0aeab", FG = "#92696c", 
                 SGP = "#cca261", GO = "#f0a325", SBAP = "#51d5a7", IWF = "#d44273", 
                 BSJP = "#424160", MSSD = "#dac370", MDCH = "#2d0cd4", CVG = "#c9edd3", 
                 SAS = "#92b1b6", CCH = "#7e22ca")
  spps.lookup <- copy(ccissr::T1)
  edatope.names <- c("Poor-subxeric", "Medium-mesic", "Rich-hygric")
  edatopes <- c("B2", "C4", "D6")
  
  # spp <- "Pl"
  # edatope <- "C4"
  # period <- "2041_2060"
  # outline <- vect("data-raw/data_tables/bc_outline.gpkg")

  if(save_png){
    if(three_panel){
      png(file=paste("./Three_Panel",spp,edatope,period,"png",sep = "."), type="cairo", units="in", width=6.5, height=2.9, pointsize=9, res=400)
    } else {
      png(file=paste("./Two_Panel",spp,edatope,period,"png",sep = "."), type="cairo", units="in", width=6.5, height=5, pointsize=12, res=300)
    }
  }
  
  par(plt=c(0,1,0,1), bg="white")
  plot(0, col="white", xaxt="n", yaxt="n", xlab="", ylab="")
  Common <- as.character(spps.lookup$EnglishName[which(spps.lookup$TreeCode==spp)])
  Latin <- as.character(spps.lookup$ScientificName[which(spps.lookup$TreeCode==spp)])
  mtext(if(spp%in%spps.lookup$TreeCode) bquote(bold(.(spp))~"-"~.(Common)) else bquote(bold(.(spp))),
        side=3, line=-2.5, adj=0.01, cex=0.9, font=2)
  if(edatope %in% edatopes) {
    mtext(paste("Site type: ", edatope, " (", edatope.names[edatope == edatopes], ")", sep=""), side=3, line=-3.5, adj=0.01, cex=0.8, font=1)
  } else {
    mtext(paste("Site type: ", edatope, sep=""), side=3, line=-3.5, adj=0.01, cex=0.8, font=1)
  }
  mtext(paste("Time period: ", period, sep=""), side=3, line=-4.5, adj=0.01, cex=0.8, font=1)
  
  ##=================================
  ###historic suitability
  dat_spp <- dbGetQuery(dbCon, sprintf("select * from cciss_res where Spp = '%s' AND FuturePeriod = '%s' AND Edatope = '%s'", spp, period, edatope)) |> as.data.table()
  dat_spp[,FeasChange := Curr - Newsuit]
  dat_spp[,Curr := as.integer(round(Curr))]
  dat_spp[is.na(Curr) | Curr > 3.5, Curr := 5]
  X <- copy(bgc_template$bgc_rast)
  values(X) <- NA

  X[dat_spp$SiteRef] <- dat_spp$Curr
  breakseq <- c(0.5,1.5,2.5,3.5,5)
  ColScheme <- c("darkgreen", "dodgerblue1", "gold2", "white")
  
  if(three_panel){
    par(plt = c(0, 0.3, 0, 0.6),new = TRUE, xpd = TRUE)
  } else {
    par(plt = c(0, 0.5, 0.05, 0.6),new = TRUE, xpd = TRUE)
  }
  
  image(X,xlab = NA,ylab = NA,bty = "n",  xaxt="n", yaxt="n",
        col=ColScheme, breaks=breakseq,asp = 1)
  plot(outline, add=T, border="black",col = NA, lwd=0.4)
  par(xpd = NA)
  legend("topleft", legend=c("1 (primary)", "2 (secondary)", "3 (tertiary)"), fill=ColScheme, bty="n", cex=0.8, title="Historical feasibility", inset=c(0,-0.3))
  # mtext(paste("(", letters[1],")", sep=""), side=3, line=-2.75, adj=0.05, cex=0.8, font=2)
  
  
  ##=================================
  ##mean feasibility change
  
  values(X) <- NA
  X2 <- copy(X)
  X3 <- copy(X)
  X[dat_spp[Curr < 3.5 | Newsuit < 3.5, SiteRef]] <- dat_spp[Curr < 3.5 | Newsuit < 3.5, FeasChange]
  X2[dat_spp[Curr > 3.5 & Newsuit < 3.5, SiteRef]] <- dat_spp[Curr > 3.5 & Newsuit < 3.5, FeasChange]
  X3[dat_spp[Curr<4 & Newsuit>3.5, SiteRef]] <- 1
  
  breakpoints <- seq(-3,3,0.5); length(breakpoints)
  labels <- c("-3","-2", "-1", "no change", "+1","+2","+3")
  ColScheme <- c("black", brewer.pal(11,"RdBu")[c(1,2,3,4)], "grey90", "grey90", brewer.pal(11,"RdBu")[c(7,8,9,10,11)]);
  ColScheme2 <- c(brewer.pal(11,"RdBu")[c(1,2,3,4,4)], "grey90", colorRampPalette(c("white", "khaki1", "gold"))(6));
  ColScheme3 <- 1
  
  if(three_panel){
    par(plt = c(0.25,0.75,0,1), xpd = TRUE, new = TRUE)
  } else {
    par(plt = c(0.25, 0.95, 0.175, 1), xpd = TRUE, new = TRUE)
  }
  
  image(X,xlab = NA,ylab = NA,bty = "n", xaxt="n", yaxt="n", col=ColScheme, breaks=breakpoints, asp = 1)
  image(X2, add=T, xlab = NA,ylab = NA,bty = "n", xaxt="n", yaxt="n", col=ColScheme2, breaks=breakpoints, asp = 1)
  image(X3, add=T, xlab = NA,ylab = NA,bty = "n", xaxt="n", yaxt="n", col=ColScheme3, asp = 1)
  plot(outline, add=T, border="black",col = NA, lwd=0.4)
  
  xl <- 1600000; yb <- 1000000; xr <- 1700000; yt <- 1700000; xadj <- 10000
  if(three_panel) {
    xl <- 325000; yb <- 900000; xr <- 400000; yt <- 1525000; xadj <- 10000
  }
  par(xpd = NA)
  y.int <- (yt-yb)/length(ColScheme)
  rect(xl+xadj,  head(seq(yb,yt,y.int),-1),  xr,  tail(seq(yb,yt,y.int),-1),  col=ColScheme)
  rect(xl-diff(c(xl+xadj, xr)),  head(seq(yb,yt,y.int),-1),  xl-xadj,  tail(seq(yb,yt,y.int),-1),  col=ColScheme2)
  rect(xl-diff(c(xl+xadj, xr)),  yb,  xl-xadj,  (yb+yt)/2,  col="white")
  text(xl-diff(c(xl+xadj, xr))/2, yb+(yt-yb)/4, "Expansion", srt=90, cex=0.85, font=1)
  text(rep(xr-10000,length(labels)),seq(yb,yt,(yt-yb)/(length(labels)-1)),labels,pos=4,cex=0.8,font=1)
  text(xl-diff(c(xl+xadj, xr))-30000, mean(c(yb,yt))-30000, paste("Mean change in feasibility", sep=""), srt=90, pos=3, cex=0.85, font=2)
  rect(xl+xadj,  yb-y.int-20000,  xr,  yb-20000,  col="black")
  text(xr, yb-y.int/2-30000, "Loss", pos=4, cex=0.8, font=1)
  

  ##=================================
  ## Summary by zone
  bgc_mapped <- as.data.frame(bgc_template$bgc_rast, cells = TRUE) 
  setDT(bgc_mapped)
  bgc_mapped[bgc_template$ids, bgc := i.bgc, on = "bgc_id"]
  bgc_mapped[, bgc_id := NULL]
  bgc_mapped[, zone := regmatches(bgc, regexpr("^[A-Z]+", bgc))]
  
  dat_spp[bgc_mapped, zone := i.zone, on = c(SiteRef = "cell")]
  
  if(three_panel){
    par(xpd=F, mar=c(4.5,2,0.1,0.1), plt = c(0.79, 0.995, 0.1, 0.275), new = TRUE, mgp=c(1.25,0.15,0))
  } else {
    par(mar=c(0,0,0,0), plt = c(0.77, 0.995, 0.001, 0.31), new = TRUE, mgp=c(1.25,0.15,0))
    plot(0, xlim=c(0,1), ylim=c(0,1), col="white", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
    
    par(mar=c(4.5,2,0.1,0.1), plt = c(0.7, 0.995, 0.08, 0.2), new = TRUE, mgp=c(1.25,0.15,0))
  }
  
  ylim=c(-3,3)
  zones_curr <- unique(dat_spp$zone)
  xlim=c(1, length(unique(dat_spp$zone)))
  z <- boxplot(FeasChange~zone, data = dat_spp, ylab="", vertical = TRUE, plot=F)
  for(i in 1:length(zones_curr)){
    temp <- dat_spp[zone == zones_curr[i], FeasChange]
    z$stats[c(1,5), i] <- quantile(temp[!is.na(temp)],c(0.05, 0.95))
  }
  bxp(z, xlim=xlim, ylim=ylim, xaxt="n", yaxt="n", xaxs="i", ylab="", pch=0,outline=FALSE)
  lines(c(0, length(unique(dat_spp$zone))+1), c(0,0), lwd=2, col="darkgrey")
  bxp(z, add=T, boxfill = zoneScheme[match(zones_curr, names(zoneScheme))], xaxt="n", yaxt="n", xaxs="i", ylab="", pch=0,outline=FALSE)
  axis(1, at=1:length(zones_curr), zones_curr, tick=F, las=2, cex.axis=0.65)
  axis(2,at=seq(ylim[1], ylim[2], 3), seq(ylim[1], ylim[2], 3), las=2, tck=0)
  mtext("Mean change in feasibility", side=3, line=0.1, adj=.975, cex=0.65, font=2)
  
  if(three_panel) {
    values(X) <- NA
    X[dat_spp[FeasChange > 0,SiteRef]] <- dat_spp[FeasChange > 0, Improve]
    X[dat_spp[FeasChange < 0,SiteRef]] <- 0 - dat_spp[FeasChange < 0, Decline]
    X[X2<0.5] <- NA # remove cells where the feasibility expansion is less than 0.5 (X2 is from the mean feasibilty panel)
    
    breakpoints <- c(seq(-100, -50,10), seq(60, 100,10));length(breakpoints)
    labels <- c("Decline", "Improve")
    ColScheme <- c(brewer.pal(11,"RdBu")[c(1:4)], "grey90", "grey90", brewer.pal(11,"RdBu")[c(8:11)]); length(ColScheme)
    
    par(plt = c(0.6, 0.95, 0.25, 1), xpd = TRUE, new = TRUE)
    image(X,xlab = NA,ylab = NA,bty = "n", xaxt="n", yaxt="n", col=ColScheme, breaks=breakpoints, asp = 1)
    plot(outline, add=T, border="black",col = NA, lwd=0.4)
    
    xl <- 1600000; yb <- 1000000; xr <- 1700000; yt <- 1700000
    rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme)
    text(rep(xr+10000,length(labels)),seq(yb,yt,(yt-yb)/(15-1))[c(3,9)],labels,pos=4,cex=0.7,font=0.7, srt=90)
    text(rep(xr-20000,length(labels)),seq(yb,yt,(yt-yb)/(15-1))[c(1,8,15)],c("100%", "50%", "100%"),pos=4,cex=0.7,font=1)
    text(xl-30000, mean(c(yb,yt))-30000, paste("Ensemble agreement\n(% of GCMs)", sep=""), srt=90, pos=3, cex=0.75, font=2)
    par(xpd=F)
  }
  
  if(save_png){
    dev.off()
  }
} 


# dat_spp[Newsuit > 3.5 & Curr <= 3, FeasChange := -10]
# dat_spp[Curr > 3.5, FeasChange := round(FeasChange) * 10]
# dat_spp[,FeasChange := round(FeasChange/0.5)*0.5]
# dat_spp[,FeasRound := round(Newsuit)]
# dat_spp[,CurrRound := round(Curr)]
# dat_spp[CurrRound > 3, CurrRound := 999]
# dat_spp[FeasRound > 3, FeasRound := 999]
# dat_spp[,AddRet := Improve]
# dat_spp[Decline > Improve, AddRet := -Decline]
# dat_spp[,AddRet := round(AddRet/20)*20]
