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
      geom_line(alpha = 0.05) +
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
      geom_line(alpha = 0.05) +
      geom_line(data = mean_spline, aes(group = ssp), linewidth = 1.4) +
      facet_wrap(~Edatopic)+
      theme_minimal() +
      labs(y="Proportion of Historic Suitable Area", x="Time Period", colour="SSPs") +
      scale_colour_manual(values = c("#E6AB02", "#D95F02", "#7F0000"))
  }
}


#' Create C.R. Mahony's bubbleplot of species expansion/persistance
#' @param persist_expand data.table. Usually created with `spp_persist_expand`
#' @param period Character. Period to create plot for. Must be one of climr::list_gcm_periods(). 
#' @param scenario Character. ssp to use for plot. Must be one of climr::list_spps(). 
#' @param species Character vector of species or "auto" to display all available species. Default "auto".
#' @param species.focal Character A single species for which to display time evolution and ensemble variation. Must be an element of the vector defined in the 'species' parameter. Default NULL.
#' @param focal.color Character A color for the focal species data points. 
#' @param edatope Character of desired edatope to use for multiple species plot, or vector of edatopes for single species.
#' @param by Character. Either "species", to display multiple species for one edatope, or "edatopes", to display multiple edatopes for one species. 
#' @param xlab  Character. x axis title.
#' @param ylab  Character. y axis title.
#' @param xlabels  logical. x axis labels.
#' @param ylabels  logical. y axis labels.
#' @param mar  numeric. plot margins.
#' @return NULL. Creates plot
#' @import data.table
#' @importFrom car dataEllipse
#' @importFrom RColorBrewer brewer.pal
#' @importFrom plotrix arctext
#' @importFrom stinepack stinterp
#' @importFrom climr list_gcm_periods
#' @export

spp_bubbleplot <- function(persist_expand, 
                           period = "2041_2060", 
                           scenario = "ssp245", 
                           species = "auto", 
                           species.focal = NULL,
                           focal.color = "lightskyblue2",
                           edatope = "C4", 
                           by = "species", 
                           xlab = "Persistence within historically suitable range",
                           ylab = "Expansion beyond historically suitable range", 
                           xlim = c(0, 1.15),
                           ylim = c(-6,4),
                           xlabels = TRUE, 
                           ylabels = TRUE,
                           mar = c(3,4,0.1,0.1)
) {
  
  if (!requireNamespace("plotrix", quietly = TRUE)) {
    stop("Package 'plotrix' is required for arctext() annotations.", call. = FALSE)
  }
  if (!requireNamespace("stinepack", quietly = TRUE)) {
    stop("Package 'stinepack' is required for stinterp() splines", call. = FALSE)
  }
  
  persist_expand <- copy(na.omit(persist_expand, col = c("Persistance","Expansion")))
  period_sel <- period
  
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
  
  period.names <- c("2001-2020", "2021-2040", "2041-2060", "2061-2080", "2081-2100")
  
  par(mar=mar, mgp=c(1.25, 0.25, 0), cex=1)
  
  plot(0, xlim=xlim, ylim=ylim, col="white", xaxt="n", yaxt="n", xlab=xlab, ylab="")
  if(xlabels) axis(1,at=seq(xlim[1], xlim[2], 0.2), labels=paste(seq(xlim[1], xlim[2], 0.2)*100,"%", sep=""), tck=0)
  if(ylabels) axis(2,at=seq(ylim[1], ylim[2]), labels=paste(round(2^(seq(ylim[1], ylim[2]))*100),"%", sep=""), las=2, tck=0)
  par(mgp=c(2.75, 0.25, 0))
  title(ylab=ylab, cex.lab=1)
  iso <- seq(0,1.2, 0.001)
  lines(1-iso, log2(iso), lty=2, lwd=2, col="darkgray")
  # if(annotate){ 
  # arctext(x = "Growing suitable range", center = c(-1, -28.7), radius = 4.6, start = 0.431*pi , cex = 0.8, stretch = 1.05, col="darkgray", font=2)
  # arctext(x = "Shrinking suitable range", center = c(-1, -29.3), radius = 4.6, start = 0.431*pi , cex = 0.8, stretch = 1.05, col="darkgray", font=2)
  # }
  #mtext(paste(edatope.names[which(edatopes==edatope)], " sites", " (", edatope, ")", sep=""), side=3, line=-1.25, adj= if(edatope=="C4") 0.025 else 0.075, cex=0.7, font=1)
  
  # log-transform Expansion values
  persist_expand[Expansion<2^(ylim[1]-1), Expansion := 2^(ylim[1]-1)] 
  persist_expand[, Expansion := log2(Expansion)]
  
  #spps <- spps[-9]
  for(i in 1:length(spps)){
    if(by == "species") {
      spp_sel <- spps[i]
      eda_sel <- edatope[1]
    } else {
      spp_sel <- species[1]
      eda_sel <- spps[i]
    }
    col.focal <- if(is.null(species.focal)) sppcolors[i] else if(spp_sel==species.focal) focal.color else "lightgray"
    col.focal2 <- if(is.null(species.focal)) "black" else if(spp_sel==species.focal) "black" else "darkgray"
    x <- persist_expand[(ssp == scenario | is.na(ssp)) & period == period_sel & Edatopic == eda_sel & spp == spp_sel, Persistance]
    y <- persist_expand[(ssp == scenario | is.na(ssp)) & period == period_sel & Edatopic == eda_sel & spp == spp_sel, Expansion]
    
    # points(x,y)
    if(length(x)>1 & var(x)>0){
      if(var(y)==0) {
        lines(range(x), range(y), col=col.focal)
      }  else {
        dataEllipse(x, y, levels=0.5, center.pch=NULL, add=T, col=col.focal, fill=T, lwd=0.5, plot.points=F)
      } 
    }
    xuse <- mean(x) 
    yuse <- mean(y)
    xuse[xuse > (xlim[2] - 0.1)] <- xlim[2] - 0.1
    yuse[yuse > (ylim[2] - 0.1)] <- ylim[2] - 0.1
    points(xuse,yuse, pch=21, bg=col.focal, cex=3, col=col.focal2)
    text(xuse,yuse, spps[i], cex=0.7, font=2, col=col.focal2)
    
  }
  
  if(!is.null(species.focal)){
    
    spp_sel <- species.focal
    
    # -----------------------
    # points for ensemble runs
    
    x <- persist_expand[(ssp == scenario | is.na(ssp)) & period == period_sel & Edatopic == eda_sel & spp == spp_sel, Persistance]
    y <- persist_expand[(ssp == scenario | is.na(ssp)) & period == period_sel & Edatopic == eda_sel & spp == spp_sel, Expansion]
    points(x,y, pch=21, bg=focal.color, cex=1)
    
    # -----------------------
    # line for ensemble mean trajectory
    
    x2 <- persist_expand[(ssp == scenario | is.na(ssp)) & Edatopic == eda_sel & spp == spp_sel, mean(Persistance), by = period][order(period), V1]
    y2 <- persist_expand[(ssp == scenario | is.na(ssp)) & Edatopic == eda_sel & spp == spp_sel, mean(Expansion), by = period][order(period), V1]
    
    # add an origin point
    x2 <- c(1, x2)
    y2 <- c(ylim[1]-1, y2)
    
    if(length(unique(sign(diff(x2))))==1 & sum(diff(x2))!=0){
      x3 <- if(unique(sign(diff(x2)))==-1) rev(x2) else x2
      y3 <- if(unique(sign(diff(x2)))==-1) rev(y2) else y2
      s <- stinterp(x3,y3, seq(min(x3),max(x3), diff(xlim)/500)) # way better than interpSpline, not prone to oscillations
      lines(s, col=1, lwd=1.5, lty=1)
    } else lines(x2, y2, col=1, lwd=1.5, lty=1)
    
    points(x2,y2, pch=21, bg=1, cex=1)
    text(x2[2],y2[2], period.names[1], pos=4, cex=0.75, font=2, offset=0.3)
    text(x2[6],y2[6], period.names[5], pos=2, cex=0.75, font=2, offset=0.3)
    
    # -----------------------
    # label
    
    j <- which(list_gcm_periods() == period_sel) + 1 #add one because there is an origin (baseline) in the vector
    points(x2[j],y2[j], pch=21, bg=focal.color, cex=3.5, col=1)
    text(x2[j],y2[j], species.focal, cex=0.8, font=2, col=1)
    
  }
  
  box()
  
}

#' Create C.R. Mahony's bubbleplot of bgc expansion/persistance
#' @param persist_expand data.table. Usually created with `bgc_persist_expand`
#' @param period Character. Period to create plot for.
#' @param scenario Character. ssp to use for plot.
#' @param xlab  Character. x axis title.
#' @param ylab  Character. y axis title.
#' @param xlabels  logical. x axis labels.
#' @param ylabels  logical. y axis labels.
#' @param mar numeric. plot margins.
#' @param plt numeric. par(plt) parameters for overlaying into an existing plot window. 
#' @return NULL. Creates plot
#' @import data.table
#' @importFrom car dataEllipse
#' @export
bgc_bubbleplot <- function(persist_expand, 
                           period = "2041_2060", 
                           scenario = "ssp245", 
                           xlab = "Climate analog within historical range",
                           ylab = "Climate analog beyond historical range", 
                           unit.focal = NULL,
                           focal.color = "lightskyblue2",
                           xlim = c(0, 1.1),
                           ylim = c(-6,4),
                           xlabels = TRUE, 
                           ylabels = TRUE,
                           mar = c(3,4,0.1,0.1), 
                           plt = NULL
) {
  
  persist_expand <- na.omit(persist_expand, col = c("Persistance","Expansion"))
  
  # colour scheme
  zone_colours <- unique(WNA_BGCs[, .(Zone, ZoneColour)])
  colScheme <- rbind(unique(WNA_BGCs[, .(BGC, SubzoneColour)]),unique(WNA_BGCs[, .(Zone, ZoneColour)]), use.names=FALSE )
  colScheme <- setNames(colScheme, c("classification", "colour"))
  
  units <- unique(persist_expand$bgc)
  period_sel <- period
  
  par(mar=mar, mgp=c(1.25, 0.25, 0), cex=1)
  period.names <- c("2001-2020", "2021-2040", "2041-2060", "2061-2080", "2081-2100")
  
  if(!is.null(plt)) par(plt = plt, new = TRUE)
  
  plot(0, xlim=xlim, ylim=ylim, col="white", xaxt="n", yaxt="n", xlab=xlab, ylab="")
  if(xlabels) axis(1,at=seq(xlim[1], xlim[2], 0.2), labels=paste(seq(xlim[1], xlim[2], 0.2)*100,"%", sep=""), tck=0)
  if(ylabels) axis(2,at=seq(ylim[1], ylim[2]), labels=paste(round(2^(seq(ylim[1], ylim[2]))*100),"%", sep=""), las=2, tck=0)
  par(mgp=c(2.75, 0.25, 0))
  title(ylab=ylab, cex.lab=1)
  iso <- seq(0,1.2, 0.001)
  lines(1-iso, log2(iso), lty=2, lwd=2, col="darkgray")
  
  persist_expand[Expansion<2^(ylim[1]-1), Expansion := 2^(ylim[1]-1)]
  persist_expand[, Expansion := log2(Expansion)]
  
  for(unit in units){
    col.focal <- if(is.null(unit.focal) || unit == unit.focal) colScheme$colour[which(colScheme$classification==unit)] else "lightgray"
    col.focal2 <- if(is.null(unit.focal) || unit == unit.focal) "black" else "darkgray"
    x <- persist_expand[(ssp == scenario | is.na(ssp)) & period == period_sel & bgc_pred == unit, Persistance]
    y <- persist_expand[(ssp == scenario | is.na(ssp)) & period == period_sel & bgc_pred == unit, Expansion]
    
    # points(x,y)
    if(length(x)>1 & var(x) > 0){
      if(var(y)==0) {
        lines(range(x, na.rm=T), range(y), col=col.focal)
      } else {
        dataEllipse(x, y, levels=0.5, center.pch=21, add=T, col=col.focal, fill=T, lwd=0.5, plot.points=F)
      } 
    }
    xuse <- mean(x) 
    yuse <- mean(y)
    xuse[xuse > (xlim[2] - 0.1)] <- xlim[2] - 0.1
    yuse[yuse > (ylim[2] - 0.1)] <- ylim[2] - 0.1
    points(xuse,yuse, pch=21, bg=col.focal, cex= 3.5, col=col.focal2)
    text(xuse,yuse, unit, cex= 0.7, font=2, col=col.focal2)
  }
  
  if(!is.null(unit.focal)){
    
    unit <- unit.focal
    
    # -----------------------
    # points for ensemble runs
    
    x <- persist_expand[(ssp == scenario) & period == period_sel & bgc_pred == unit, Persistance]
    y <- persist_expand[(ssp == scenario) & period == period_sel & bgc_pred == unit, Expansion]
    points(x,y, pch=21, bg=focal.color, cex=1)
    
    # -----------------------
    # line for ensemble mean trajectory
    
    x2 <- persist_expand[(ssp == scenario) & bgc_pred == unit, mean(Persistance), by = period][order(period), V1]
    y2 <- persist_expand[(ssp == scenario) & bgc_pred == unit, mean(Expansion), by = period][order(period), V1]
    
    # add an origin point
    x2 <- c(1, x2)
    y2 <- c(ylim[1]-1, y2)
    
    if(length(unique(sign(diff(x2))))==1 & sum(diff(x2))!=0){
      x3 <- if(unique(sign(diff(x2)))==-1) rev(x2) else x2
      y3 <- if(unique(sign(diff(x2)))==-1) rev(y2) else y2
      s <- stinterp(x3,y3, seq(min(x3),max(x3), diff(xlim)/500))
      lines(s, col=1, lwd=1.5, lty=1)
    } else lines(x2, y2, col=1, lwd=1.5, lty=1)
    
    points(x2,y2, pch=21, bg=1, cex=1)
    text(x2[2],y2[2], period.names[1], pos=4, cex=0.75, font=2, offset=0.3)
    text(x2[6],y2[6], period.names[5], pos=2, cex=0.75, font=2, offset=0.3)
    
    # -----------------------
    # label
    
    j <- which(list_gcm_periods() == period_sel) + 1 #add one because there is an origin (baseline) in the vector
    points(x2[j],y2[j], pch=21, bg=focal.color, cex=3.5, col=1)
    text(x2[j],y2[j], unit, cex=0.8, font=2, col=1)
    
  }
  box()
}

#' Create alluvial/stacked bar plot of projected area by zone/subzone
#' @param spp Character. Single species to use for plot
#' @param edatope Character. Single edatopic position for plot (e.g., "C4")
#' @param bgc_template List containing SpatRaster of BGCs and id table, or data.table (must have columns `cell`,`BGC`). Usually created using `make_bgc_template`
#' @param fractional Logical. Use fractional (suitability based) values for calculations?
#' @param save_png Logical. Save plot to png? Default `TRUE`. If `FALSE`, creates plot on default plotting device.
#' @param by_zone Logical. Plot by zone or subzone? Defaults to `TRUE`
#' @param base_folder Base folder to read results from.
#' @return NULL. Creates plot
#' @import data.table duckdb ggplot2 scales
#' @importFrom ggalluvial geom_alluvium
#' @export

plot_spparea <- function(dbCon, 
                         spp, 
                         edatope, 
                         fractional, 
                         by_zone = TRUE, 
                         save_png = TRUE, 
                         width = 7, 
                         height = 5.5, 
                         res = 300) {
  cciss_spp <- dbGetQuery(dbCon, sprintf("select * from cciss_res where Spp = '%s' AND Edatope = '%s'", spp, edatope)) |> as.data.table()
  bgc_mapped <- dbGetQuery(dbCon, "select * from bgc_points") |> as.data.table()
  
  if (save_png) {
    png(file = paste("./alluv", spp, edatope, "png", sep = "."), type = "cairo", units = "in", width = width, height = height, res = res)
    on.exit(dev.off(), add = TRUE)
  }
  
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
  
  yrs <- sort(unique(as.integer(cciss_sum$Year)))
  idx <- seq_along(yrs)
  
  midpts <- head(idx, -1) + 0.5
  
  bars <- data.frame(
    xmin = midpts - 0.3,
    xmax = midpts + 0.3,
    ymin = -Inf,
    ymax = Inf
  )
  
  # colour scheme
  zone_colours <- unique(WNA_BGCs[, .(Zone, ZoneColour)])
  colScheme <- rbind(unique(WNA_BGCs[, .(BGC, SubzoneColour)]),unique(WNA_BGCs[, .(Zone, ZoneColour)]), use.names=FALSE )
  colScheme <- setNames(colScheme, c("classification", "colour"))
  col_final <- colScheme$colour
  names(col_final) <- colScheme$classification
  
  spat_res <- dbGetQuery(dbCon, "select * from spatial_res") |> as.data.table()
  if(spat_res$projected[1]){
    cellarea <- (spat_res$y_res/1000) * (spat_res$x_res/1000)
  } else {
    cellarea <- (spat_res$y_res*111)*(spat_res$x_res*111*cos(mean(c(spat_res$y_min,spat_res$y_max)) * pi / 180))
    warning("Input data is not in projected crs. Cell area calculations will be approximate.")
  }
  
  cciss_sum[, SppArea := SppArea * cellarea]
  
  year_levels <- sort(as.character(cciss_sum$Year))
  
  # order zones by change (last - first): most decline at bottom
  delta_by_zone <- cciss_sum[SppArea > 0,
                             .(SppArea = sum(SppArea, na.rm = TRUE)),
                             by = .(zone, Year)]
  delta_by_zone <- dcast(delta_by_zone, zone ~ Year, value.var = "SppArea", fill = 0)
  delta_by_zone[, delta := get(tail(year_levels, 1)) - get(head(year_levels, 1))]
  zone_order <- delta_by_zone[order(delta), zone]
  zone_order <- rev(zone_order) # should put declines at bottom, increases at top
  cciss_sum[, zone := factor(zone, levels = zone_order)]
  cciss_sum <- cciss_sum[!is.na(zone) & SppArea > 0]
  
  ## fill in with zeros
  grid <- CJ(
    zone   = sort(unique(cciss_sum$zone)),
    Year = sort(unique(cciss_sum$Year)),
    unique = TRUE
  )
  
  # left join onto grid, then fill missing with 0
  cciss_sum_full <- cciss_sum[grid, on = .(zone, Year)]
  cciss_sum_full[is.na(SppArea), SppArea := 0]
  
  # Plot
  ggplot(cciss_sum_full, aes(x = Year, y = SppArea, fill = zone)) +
    geom_alluvium(aes(alluvium = zone), alpha= 1, color = "black") +
    
    geom_rect(
      data = bars,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "white", alpha = 0.8
    ) +
    
    #add lines bordering each column
    geom_col(
      data = cciss_sum_full |>
        dplyr::group_by(Year) |>
        dplyr::summarise(total = sum(SppArea), .groups = "drop"),
      aes(x = Year, y = total),
      position = "stack",
      width = 0.4, 
      fill = NA,
      color = "black",
      linewidth = 1.5,
      inherit.aes = FALSE
    ) +
    
    theme_classic(base_size = 12) +
    
    #legend matches stack order + drop unplotted zones
    scale_fill_manual(
      values = col_final,
      breaks = zone_order,
      limits = zone_order,
      drop = TRUE
    ) +
    
    scale_x_discrete(labels=c("1961" = "1961-90", "2001" = "2001-20", "2021" = "2021-40",
                              "2041" = "2041-60", "2061" = "2061-80", "2081" = "2081-2100")) +
    scale_y_continuous(labels = scales::comma, expand=c(0,0)) +
    labs(y="Environmentally Suitable Area (Km^2)",x="Time period", fill = "BGC zone") +
    theme(axis.ticks.x = element_blank())
  
}

#' Create alluvial/stacked bar plot of projected area by zone/subzone
#' @param dat data.table of species area data from 'alluvial_area' table
#' @param spp Character. Single species to use for plot
#' @param edatope Character. Single edatopic position for plot (e.g., "C4")
#' @param by_zone Logical. Plot by zone or subzone? Defaults to `TRUE`
#' @return NULL. Creates plot
#' @import data.table ggplot2 scales
#' @importFrom ggalluvial geom_alluvium
#' @export
plot_alluvial <- function(dat, spp, edatope, by_zone = T, cellarea = 4) {
  dat <- dat[Spp == spp & Edatope == edatope,]
  dat[, Year := as.factor(substr(FuturePeriod,1,4))]
  
  yrs <- sort(unique(as.integer(dat$Year)))
  idx <- seq_along(yrs)
  
  midpts <- head(idx, -1) + 0.5
  
  bars <- data.frame(
    xmin = midpts - 0.3,
    xmax = midpts + 0.3,
    ymin = -Inf,
    ymax = Inf
  )
  
  # colour scheme
  zone_colours <- unique(WNA_BGCs[, .(Zone, ZoneColour)])
  colScheme <- rbind(unique(WNA_BGCs[, .(BGC, SubzoneColour)]),unique(WNA_BGCs[, .(Zone, ZoneColour)]), use.names=FALSE )
  colScheme <- setNames(colScheme, c("classification", "colour"))
  colScheme <- colScheme[names(colScheme) %in% unique(dat[SppArea > 0,bgc])]
  
  dat[, SppArea := SppArea * cellarea]
  
  year_levels <- sort(as.character(dat$Year))
  
  # order zones by change (last - first): most decline at bottom
  delta_by_zone <- dat[,.(SppArea = sum(SppArea, na.rm = TRUE)),
                             by = .(bgc, Year)]
  delta_by_zone <- dcast(delta_by_zone, bgc ~ Year, value.var = "SppArea", fill = 0)
  delta_by_zone[, delta := get(tail(year_levels, 1)) - get(head(year_levels, 1))]
  zone_order <- delta_by_zone[order(delta), bgc]
  zone_order <- rev(zone_order) # should put declines at bottom, increases at top
  zone_order <- zone_order[zone_order %in% unique(dat[SppArea > 0,bgc])]
  dat[, bgc := factor(bgc, levels = zone_order)]
  dat <- dat[!is.na(bgc) & SppArea > 0]
  
  ## fill in with zeros
  grid <- CJ(
    bgc   = sort(unique(dat$bgc)),
    Year = sort(unique(dat$Year)),
    unique = TRUE
  )
  
  # left join onto grid, then fill missing with 0
  cciss_sum_full <- dat[grid, on = .(bgc, Year)]
  cciss_sum_full[is.na(SppArea), SppArea := 0]
  cciss_sum_full[,bgc := factor(as.character(bgc), levels = zone_order)]
  
  # Plot
  ggplot(cciss_sum_full, aes(x = Year, y = SppArea, fill = bgc)) +
    geom_alluvium(aes(alluvium = bgc), alpha= 1, color = "black") +
    
    geom_rect(
      data = bars,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "white", alpha = 0.8
    ) +
    
    #add lines bordering each column
    geom_col(
      data = cciss_sum_full |>
        dplyr::group_by(Year) |>
        dplyr::summarise(total = sum(SppArea), .groups = "drop"),
      aes(x = Year, y = total),
      position = "stack",
      width = 0.4, 
      fill = NA,
      color = "black",
      linewidth = 1.5,
      inherit.aes = FALSE
    ) +
    
    theme_classic(base_size = 12) +
    
    #legend matches stack order + drop unplotted zones
    scale_fill_manual(
      values = colScheme,
      breaks = zone_order,
      limits = zone_order,
      drop = TRUE
    ) +
    
    scale_x_discrete(labels=c("1961" = "1961-90", "2001" = "2001-20", "2021" = "2021-40",
                              "2041" = "2041-60", "2061" = "2061-80", "2081" = "2081-2100")) +
    scale_y_continuous(labels = scales::comma, expand=c(0,0)) +
    labs(y=expression("Environmentally suitable area (km"^2*")"), x="Time period", fill = "BGC zone") +
    theme(axis.ticks.x = element_blank())
  
}

#' Plot map of reference suitability for given species and edatope
#' @description
#' Note that currently, this function only works correctly if the preceeding analysis has been done using an Albers grid.
#' @param con duckdb database connection
#' @param raster_template SpatRaster to use as template for plotting.
#' @param outline SpatVector of aoi boundary.
#' @param spp Character. Species code to create plot for
#' @param edatope Character. Edatope code to create plot for (e.g., "C4")
#' @param save_png Logical. Save plot to png? Default `FALSE`. If `FALSE`, creates plot on default plotting device.
#' @return NULL. Creates plot
#' @import data.table duckdb terra
#' @export
plot_reference_suitability <- function(con,
                                       raster_template,
                                       outline,
                                       spp, 
                                       edatope = "C4",
                                       save_png = FALSE) {
  if(save_png){
    png(file=paste("./ReferenceSuitability",spp,edatope,"png",sep = "."), type="cairo", units="in", width=4, height=5, pointsize=12, res=300)
  }

  ##=================================
  ###historic suitability
  dat_spp <- dbGetQuery(con, sprintf("select * from cciss_res where Spp = '%s' AND FuturePeriod = '2041_2060' AND Edatope = '%s'", spp, edatope)) |> as.data.table()
  dat_spp[,Curr := as.integer(round(Curr))]
  dat_spp[is.na(Curr) | Curr > 3.5, Curr := 4]
  dat_spp[is.na(Newsuit) | Newsuit > 3.5, Newsuit := 4]
  dat_spp[,FeasChange := Curr - Newsuit]
  X <- copy(raster_template)
  values(X) <- NA
  
  X[dat_spp$SiteRef] <- dat_spp$Curr
  breakseq <- c(0.5,1.5,2.5,3.5,5)
  colScheme <- c("darkgreen", "dodgerblue1", "gold2", "white")
  
  image(X,xlab = NA,ylab = NA,bty = "n",  xaxt="n", yaxt="n",
        col=colScheme, breaks=breakseq,asp = 1)
  terra::plot(outline, add=T, border="black",col = NA, lwd=0.4)
  par(xpd = NA)
  legend("topleft", legend = c("E1 (high)", "E2 (moderate)", "E3 (low)"), fill=colScheme, bty="n", cex=0.8, title="Historical suitability", inset=c(0,-0.3))
  
  if(save_png){
    dev.off()
  }
}

#' Create C.R. Mahony's two-panel map plot of historic suitability, suitability change, and a boxplot of change by zone
#' @description
#' Note that currently, this function only works correctly if the preceeding analysis has been done using an Albers grid.
#' @param dbCon duckdb database connection
#' @param bgc_template List containing SpatRaster of BGCs and id table. Usually created using `make_bgc_template`
#' @param outline SpatVector of aoi boundary.
#' @param spp Character. Species code to create plot for
#' @param edatope Character. Edatope code to create plot for (e.g., "C4")
#' @param period Character. Period to create plot for.
#' @param save_png Logical. Save plot to png? Default `TRUE`. If `FALSE`, creates plot on default plotting device.
#' @param panel_labels Logical. Add a manuscript-style label to each panel of the plot.
#' @return NULL. Creates plot
#' @import data.table duckdb terra
#' @export
plot_SuitabilityChangeMap <- function(dbCon, 
                                      bgc_template, 
                                      outline, 
                                      spp = "Fd", 
                                      edatope = "C4", 
                                      period = "2041_2060", 
                                      three_panel = FALSE, 
                                      save_png = TRUE,
                                      panel_labels = TRUE
) 
{
  
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
  dat_spp[,Curr := as.integer(round(Curr))]
  dat_spp[is.na(Curr) | Curr > 3.5, Curr := 4]
  dat_spp[is.na(Newsuit) | Newsuit > 3.5, Newsuit := 4]
  dat_spp[,FeasChange := Curr - Newsuit]
  X <- copy(bgc_template$bgc_rast)
  values(X) <- NA
  
  X[dat_spp$SiteRef] <- dat_spp$Curr
  breakseq <- c(0.5,1.5,2.5,3.5,5)
  colScheme <- c("darkgreen", "dodgerblue1", "gold2", "white")
  
  if(three_panel){
    par(plt = c(0, 0.3, 0, 0.6),new = TRUE, xpd = TRUE)
  } else {
    par(plt = c(0, 0.5, 0.05, 0.6),new = TRUE, xpd = TRUE)
  }
  
  image(X,xlab = NA,ylab = NA,bty = "n",  xaxt="n", yaxt="n",
        col=colScheme, breaks=breakseq,asp = 1)
  terra::plot(outline, add=T, border="black",col = NA, lwd=0.4)
  par(xpd = NA)
  legend("topleft", legend = c("E1 (high)", "E2 (moderate)", "E3 (low)"), fill=colScheme, bty="n", cex=0.8, title="Historical suitability", inset=c(0,-0.3))
  
  if(panel_labels) mtext(paste("(", letters[1],")", sep=""), side=3, line=-8, adj=0.05, cex=0.8, font=2)
  
  
  ##=================================
  ##mean feasibility change
  
  values(X) <- NA
  X2 <- copy(X)
  X3 <- copy(X)
  X[dat_spp[Curr < 3.5 | Newsuit < 3.5, SiteRef]] <- dat_spp[Curr < 3.5 | Newsuit < 3.5, FeasChange]
  X2[dat_spp[Curr > 3.5 & Newsuit < 3.5, SiteRef]] <- dat_spp[Curr > 3.5 & Newsuit < 3.5, FeasChange]
  X3[dat_spp[Curr<4 & Newsuit>3.5, SiteRef]] <- 1
  
  breakpoints <- seq(-3,3,0.5); length(breakpoints)
  labels <- c("-3","-2", "-1", "0", "+1","+2","+3")
  colScheme <- c("black", brewer.pal(11,"RdBu")[c(1,2,3,4)], "grey90", "grey90", brewer.pal(11,"RdBu")[c(7,8,9,10,11)]);
  colScheme2 <- c(brewer.pal(11,"RdBu")[c(1,2,3,4,4)], "grey90", colorRampPalette(c("white", "khaki1", "gold"))(6));
  colScheme3 <- 1
  
  if(three_panel){
    par(plt = c(0.25,0.75,0,1), xpd = TRUE, new = TRUE)
  } else {
    par(plt = c(0.25, 0.95, 0.175, 1), xpd = TRUE, new = TRUE)
  }
  
  image(X,xlab = NA,ylab = NA,bty = "n", xaxt="n", yaxt="n", col=colScheme, breaks=breakpoints, asp = 1)
  image(X2, add=T, xlab = NA,ylab = NA,bty = "n", xaxt="n", yaxt="n", col=colScheme2, breaks=breakpoints, asp = 1)
  image(X3, add=T, xlab = NA,ylab = NA,bty = "n", xaxt="n", yaxt="n", col=colScheme3, asp = 1)
  terra::plot(outline, add=T, border="black",col = NA, lwd=0.4)
  
  xl <- 1600000; yb <- 1000000; xr <- 1700000; yt <- 1700000; xadj <- 10000
  if(three_panel) {
    xl <- 325000; yb <- 900000; xr <- 400000; yt <- 1525000; xadj <- 10000
  }
  par(xpd = NA)
  y.int <- (yt-yb)/length(colScheme)
  rect(xl+xadj,  head(seq(yb,yt,y.int),-1),  xr,  tail(seq(yb,yt,y.int),-1),  col=colScheme)
  rect(xl-diff(c(xl+xadj, xr)),  head(seq(yb,yt,y.int),-1),  xl-xadj,  tail(seq(yb,yt,y.int),-1),  col=colScheme2)
  rect(xl-diff(c(xl+xadj, xr)),  yb,  xl-xadj,  (yb+yt)/2,  col="white")
  text(xl-diff(c(xl+xadj, xr))/2, yb+(yt-yb)/4, "Expansion", srt=90, cex=0.85, font=1)
  text(rep(xr-10000,length(labels)),seq(yb,yt,(yt-yb)/(length(labels)-1)),labels,pos=4,cex=0.8,font=1)
  text(xl-diff(c(xl+xadj, xr))-30000, mean(c(yb,yt))-30000, paste("Mean change in suitability", sep=""), srt=90, pos=3, cex=0.85, font=2)
  rect(xl+xadj,  yb-y.int-20000,  xr,  yb-20000,  col="black")
  text(xr, yb-y.int/2-30000, "Loss", pos=4, cex=0.8, font=1)
  
  if(panel_labels) mtext(paste("(", letters[2],")", sep=""), side=3, line=1, adj=0.08, cex=0.8, font=2)
  
  ##=================================
  ## Summary by zone
  
  zone_order <- c("CDF", "CWH", "MH", "ESSF", "MS", "IDF", "PP", "BG", "ICH", "SBPS", "SBS", "BWBS", "SWB", "CMA", "IMA", "BAFA")
  
  bgc_mapped <- as.data.frame(bgc_template$bgc_rast, cells = TRUE) 
  setDT(bgc_mapped)
  bgc_mapped[bgc_template$ids, bgc := i.bgc, on = "bgc_id"]
  bgc_mapped[, bgc_id := NULL]
  bgc_mapped[, zone := factor(regmatches(bgc, regexpr("^[A-Z]+", bgc)), levels = zone_order)]
  
  dat_spp[bgc_mapped, zone := i.zone, on = c(SiteRef = "cell")]
  
  if(three_panel){
    par(xpd=F, mar=c(4.5,2,0.1,0.1), plt = c(0.79, 0.995, 0.1, 0.275), new = TRUE, mgp=c(1.25,0.15,0))
  } else {
    par(mar=c(0,0,0,0), plt = c(0.77, 0.995, 0.001, 0.31), new = TRUE, mgp=c(1.25,0.15,0))
    plot(0, xlim=c(0,1), ylim=c(0,1), col="white", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
    
    par(mar=c(4.5,2,0.1,0.1), plt = c(0.7, 0.995, 0.08, 0.2), new = TRUE, mgp=c(1.25,0.15,0))
  }
  
  ylim=c(-3,3)
  zones_curr <- levels(dat_spp$zone)[levels(dat_spp$zone) %in% dat_spp$zone]
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
  mtext("Mean change in suitability", side=3, line=0.1, adj=.975, cex=0.65, font=2)
  
  if(panel_labels) mtext(paste("(", letters[3],")", sep=""), side=3, line=1, adj=0.975, cex=0.8, font=2)
  
  if(three_panel) {
    values(X) <- NA
    X[dat_spp[FeasChange > 0,SiteRef]] <- dat_spp[FeasChange > 0, Improve]
    X[dat_spp[FeasChange < 0,SiteRef]] <- 0 - dat_spp[FeasChange < 0, Decline]
    X[X2<0.5] <- NA # remove cells where the feasibility expansion is less than 0.5 (X2 is from the mean feasibilty panel)
    
    breakpoints <- c(seq(-100, -50,10), seq(60, 100,10));length(breakpoints)
    labels <- c("Decline", "Improve")
    colScheme <- c(brewer.pal(11,"RdBu")[c(1:4)], "grey90", "grey90", brewer.pal(11,"RdBu")[c(8:11)]); length(colScheme)
    
    par(plt = c(0.6, 0.95, 0.25, 1), xpd = TRUE, new = TRUE)
    image(X,xlab = NA,ylab = NA,bty = "n", xaxt="n", yaxt="n", col=colScheme, breaks=breakpoints, asp = 1)
    terra::plot(outline, add=T, border="black",col = NA, lwd=0.4)
    
    xl <- 1600000; yb <- 1000000; xr <- 1700000; yt <- 1700000
    rect(xl,  head(seq(yb,yt,(yt-yb)/length(colScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(colScheme)),-1),  col=colScheme)
    text(rep(xr+10000,length(labels)),seq(yb,yt,(yt-yb)/(15-1))[c(3,9)],labels,pos=4,cex=0.7,font=0.7, srt=90)
    text(rep(xr-20000,length(labels)),seq(yb,yt,(yt-yb)/(15-1))[c(1,8,15)],c("100%", "50%", "100%"),pos=4,cex=0.7,font=1)
    text(xl-30000, mean(c(yb,yt))-30000, paste("Ensemble agreement\n(% of GCMs)", sep=""), srt=90, pos=3, cex=0.75, font=2)
    par(xpd=F)
    
    if(panel_labels) mtext(paste("(", letters[4],")", sep=""), side=3, line=-3.25, adj=0.1, cex=0.8, font=2)
    
  }
  
  if(save_png){
    dev.off()
  }
} 

#' Save raster of predicted suitability for use in e.g. CCISS Spatial
#' @param con duckdb connection
#' @param spp Species to plot
#' @param period Character. Period to plot for.
#' @param edatope Character. Edatope position (e.g. C4)
#' @param raster_template SpatRaster to use as template for plotting. 
#' @param save_location Character. Path to save raster in. Will be created if it doesn't exist. 
#' @param raw Logical. Save raw values, or colourised rgb rasters? Default `FALSE`
#' @param plot Logical. Create plot as well as saving raster? Default `FALSE`
#' @import data.table
#' @importFrom glue glue_sql
#' @importFrom DBI dbGetQuery
#' @importFrom terra  colorize plotRGB writeRaster as.int
#' @export
plot_predicted_suitability <- function(con, spp, period, edatope, raster_template, save_location = "suit_rasters", raw = FALSE, plot = FALSE) {
  if(!dir.exists(save_location)) dir.create(save_location)
  dat_spp <- dbGetQuery(con, glue_sql("select * from cciss_res where Spp = {spp} AND FuturePeriod = {period} AND Edatope = {edatope}", .con = con)) |> as.data.table()
  
  breakpoints.suit <- c(1,2,3,999)
  palette.suit <-   c("#006400", "#1E90FF", "#EEC900", "#FFFFFF")
  suit_cols <- data.table(value = breakpoints.suit,Colour = palette.suit)
  dat_spp[,FeasRound := round(Newsuit)]
  dat_spp[FeasRound > 3, FeasRound := 999]
  
  final_dem <- copy(raster_template)
  values(final_dem) <- NA
  
  if(!raw) {
    final_dem[!is.na(raster_template)] <- 999
    final_dem[dat_spp$SiteRef] <- dat_spp$FeasRound
    coltab(final_dem) <- suit_cols
    final_rgb <- colorize(final_dem, to = "rgb", alpha = TRUE)
    if(plot) plotRGB(final_rgb)
    writeRaster(final_rgb, paste0(save_location,"/NewFeas_",period,"_",edatope,"_",spp,".tif"), overwrite = T)
  } else {
    final_dem[dat_spp$SiteRef] <- dat_spp$FeasRound
    final_dem[final_dem == 999] <- NA
    final_dem <- as.int(final_dem * 10)
    writeRaster(final_dem, paste0(save_location,"/FeasibilityRaw_",period,"_",edatope,"_",spp,".tif"),overwrite = T, datatype = "INT2U")
  }
  
}

#' Save raster of predicted suitability change for use in e.g. CCISS Spatial
#' @param con duckdb connection
#' @param spp Species to plot
#' @param period Character. Period to plot for.
#' @param edatope Character. Edatope position (e.g. C4)
#' @param raster_template SpatRaster to use as template for plotting. 
#' @param save_location Character. Path to save raster in. Will be created if it doesn't exist. 
#' @param raw Logical. Save raw values, or colourised rgb rasters? Default `FALSE`
#' @param plot Logical. Create plot as well as saving raster? Default `FALSE`
#' @import data.table
#' @importFrom glue glue_sql
#' @importFrom DBI dbGetQuery
#' @importFrom terra  colorize plotRGB writeRaster as.int
#' @export
plot_suitability_change <- function(con, spp, period, edatope, raster_template, save_location = "suit_rasters", raw = FALSE, plot = FALSE) {
  if(!dir.exists(save_location)) dir.create(save_location)
  dat_spp <- dbGetQuery(con, glue_sql("select * from cciss_res where Spp = {spp} AND FuturePeriod = {period} AND Edatope = {edatope}", .con = con)) |> as.data.table()
  
  breakpoints.change <- c(c(seq(-2.5,2.5,0.5),-10,10,20,30) + 15, 999)
  palette.change <- c(brewer.pal(11,"RdBu")[c(1,2,3,4,5,6)], brewer.pal(11,"RdBu")[c(7,8,9,10,11)],"#000000", brewer.pal(9,"YlOrRd")[1:3],"#FFFFFF") # nolint
  change_cols <- data.table(value = breakpoints.change, Colour = palette.change)
  change_cols[value == 15, Colour := "#DFDFDF"]
  
  dat_spp <- dat_spp[Curr < 3.5 | Newsuit < 3.5,]
  dat_spp[,FeasChange := Curr - Newsuit]
  dat_spp[Newsuit > 3.5 & Curr <= 3, FeasChange := -10]
  dat_spp[Curr > 3.5, FeasChange := round(FeasChange) * 10]
  dat_spp[,FeasChange := round(FeasChange/0.5)*0.5]
  dat_spp[,FeasRound := round(Newsuit)]
  dat_spp[,CurrRound := round(Curr)]
  dat_spp[CurrRound > 3, CurrRound := 999]
  dat_spp[FeasRound > 3, FeasRound := 999]

  final_dem <- copy(raster_template)
  values(final_dem) <- NA
  
  if(!raw) {
    values(final_dem) <- NA
    final_dem[!is.na(raster_template)] <- 999
    final_dem[dat_spp$SiteRef] <- dat_spp$FeasChange + 15
    final_rgb <- subst(final_dem, change_cols$value, t(col2rgb(change_cols$Colour,alpha = TRUE)),names = c("red","green", "blue","alpha"))
    if(plot) plotRGB(final_rgb)
    writeRaster(final_rgb, paste0(save_location,"/MeanChange_",period,"_",edatope,"_",spp,".tif"), overwrite = T)
  } else {
    final_dem[dat_spp$SiteRef] <- dat_spp$FeasChange
    final_dem[final_dem == 999] <- NA
    final_dem <- as.int(final_dem * 10)
    writeRaster(final_dem, paste0(save_location,"/MeanChangeRaw_",period,"_",edatope,"_",spp,".tif"),overwrite = T, datatype = "INT4S")  
  }
}


#' Save raster of predicted BGC map
#' @param con duckdb connection
#' @param period Character. Period to plot for.
#' @param plot_ensemble Logical. Plot vote winner over ensemble of runs?
#' @param gcm,ssp,run Character. Only used if `plot_ensemble = FALSE`
#' @param raster_template SpatRaster to use as template for plotting. 
#' @param by_zone Logical. Plot by zone or subzone variant? Default `FALSE`
#' @param save_location Character. Path to save raster in. Will be created if it doesn't exist. 
#' @param raw Logical. Save raw values, or colourised rgb rasters? Default `FALSE`
#' @param plot Logical. Create plot as well as saving raster? Default `FALSE`
#' @import data.table
#' @importFrom glue glue_sql
#' @importFrom DBI dbGetQuery
#' @importFrom terra  colorize plotRGB writeRaster as.int
#' @export
plot_bgc <- function(con, period, plot_ensemble, raster_template, gcm = NULL, ssp = NULL, run = NULL,
                      by_zone = FALSE, 
                     save_location = "bgc_rasters", raw = FALSE, plot = FALSE) {
  if(!dir.exists(save_location)) dir.create(save_location)
  
  if(plot_ensemble){
    if(by_zone) table_name <- DBI::SQL("ensemble_preds_zone")
    else table_name <- DBI::SQL("ensemble_preds")
    dat <- dbGetQuery(con, glue_sql("select cellnum, bgc_pred 
                                    from {table_name} 
                                    where period = {period}", .con = con)) |> as.data.table()
  } else {
    if(by_zone) sel_bgc <- DBI::SQL("regexp_extract(bgc_pred, '^[A-Z]+')")
    else sel_bgc <- DBI::SQL("bgc_pred")
    if(period == "1961_1990") {
      dat <- dbGetQuery(con, glue_sql("select cellnum, {sel_bgc} AS bgc_pred 
                                    from bgc_raw 
                                    where period = {period}", .con = con)) |> as.data.table()
    } else {
      dat <- dbGetQuery(con, glue_sql("select cellnum, {sel_bgc} AS bgc_pred 
                                    from bgc_raw 
                                    where period = {period}", .con = con)) |> as.data.table()
    } else {
      dat <- dbGetQuery(con, glue_sql("select cellnum, {sel_bgc} AS bgc_pred 
                                    from bgc_raw_runs 
                                    where ssp = {ssp} and 
                                    gcm = {gcm} and 
                                    run = {run} and 
                                    period = {period}", .con = con)) |> as.data.table()
    }
    
    }
  }
  if(by_zone) {
    col_use <- WNA_BGCs[,.(Zone, ZoneColour)]
  } else {
    col_use <- WNA_BGCs[,.(BGC, SubzoneColour)]
  }
  col_use <- unique(col_use)
  setnames(col_use, c("classification", "colour"))
  
  final_dem <- copy(raster_template)
  values(final_dem) <- NA
  if(plot_ensemble) gcm_ids <- "Ensemble"
  else gcm_ids <- paste(gcm, ssp, run, sep = "_")
  zonesz <- "Subzone"
  if(by_zone) zonesz <- "Zone"
  
  if(!raw) {
    dat[,bgc_id := as.numeric(as.factor(bgc_pred))]
    final_dem[dat$cellnum] <- dat$bgc_id
    bgc_id <- unique(dat[,.(bgc_pred,bgc_id)])
    bgc_id[col_use, colour := i.colour, on = c(bgc_pred = "classification")]
    coltab(final_dem) <- bgc_id[,.(bgc_id,colour)]
    rgbbgc <- colorize(final_dem, to = "rgb", alpha = T)
    if(plot) plotRGB(rgbbgc)
    writeRaster(rgbbgc, paste0(save_location,"/bgc_",gcm_ids,"_",period,"_",zonesz, ".tif"), overwrite=TRUE)  
  } else {
    if(by_zone) bgc_id <- data.table(bgc_pred = sort(WNA_BGCs$Zone)) |> unique()
    else bgc_id <- data.table(bgc_pred = sort(WNA_BGCs$BGC)) |> unique()
    bgc_id[,bgc_id := seq_along(bgc_pred)]
    dat[bgc_id, bgc_id := i.bgc_id, on = "bgc_pred"]
    final_dem[dat$cellnum] <- dat$bgc_id
    writeRaster(final_dem, paste0(save_location,"/bgcRaw_",gcm_ids,"_",period,"_",zonesz, ".tif"), overwrite=TRUE, datatype = "INT2U")  
  }
}


#' Save raster of actual BGC map
#' @param raster_template List returned from `create_bgc_template()`
#' @param by_zone Logical. Plot by zone or subzone variant? Default `FALSE`
#' @param save_location Character. Path to save raster in. Will be created if it doesn't exist. 
#' @param plot Logical. Create plot as well as saving raster? Default `FALSE`
#' @import data.table
#' @importFrom glue glue_sql
#' @importFrom DBI dbGetQuery
#' @importFrom terra  colorize plotRGB writeRaster as.int
#' @export
plot_bgc_mapped <- function(raster_template, by_zone = FALSE, 
                            save_location = "bgc_rasters", plot = FALSE) {
  if(!dir.exists(save_location)) dir.create(save_location)
  
  bgc_ids <- raster_template$ids
  if(by_zone) {
    col_use <- WNA_BGCs[,.(Zone, ZoneColour)]
    bgc_ids[,bgc := sub("^([A-Z]+).*", "\\1", bgc)]
  } else {
    col_use <- WNA_BGCs[,.(BGC, SubzoneColour)]
  }
  col_use <- unique(col_use)
  setnames(col_use, c("bgc", "colour"))
  bgc_ids[col_use, colour := i.colour, on = "bgc"]
  
  final_dem <- copy(raster_template$bgc_rast)
  
  zonesz <- "Subzone"
  if(by_zone) zonesz <- "Zone"
  
  coltab(final_dem) <- bgc_ids[,.(bgc_id,colour)]
  rgbbgc <- colorize(final_dem, to = "rgb", alpha = T)
  if(plot) plotRGB(rgbbgc)
  writeRaster(rgbbgc, paste0(save_location,"/bgc_mapped_", zonesz, ".tif"), overwrite=TRUE)  
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


#' Create a map of biogeoclimatic units with the standard BGC colour scheme
#' @description
#' This is a utility function that plots a basic map
#' @param X SpatRaster. A template raster.
#' @param dat data.table with first column indicating template raster cell ID and second column the BGC subzone/variant. 
#' @param zone logical. Plot with BGC zone colours. 
#' @param boundary SpatVector. Optional study area boundary.
#' @param mask SpatVector. Optional spatial mask
#' @param legend logical. Plot a legend. 
#' @param mask_alpine logical. Plot alpine BGC units as white.
#' @param label_exotic numeric. Minimum number of grid cells required to plot a label for exotic (non-BC) BGC units. 
#' @param q_exotic numeric. Quantile of cell numbers for positioning the label for exotic BGC units. Acceptable values are between 0 and 1. 
#' @param title character. Title for the plot, left justified at top. 
#' @param add logical. Add the map to an existing plot.
#' @import data.table terra
#' @export

bgc_map <- function(X, 
                    dat, 
                    zone = TRUE, 
                    boundary = NULL, 
                    mask = NULL, 
                    legend = FALSE, 
                    mask_alpine = TRUE, 
                    mask_ocean = TRUE, 
                    label_exotic = NULL, 
                    q_exotic = 0.5, 
                    title = "",
                    add = FALSE 
){
  
  if(mask_alpine){
    zones.bc <- c("BG", "BWBS", "CDF", "CWH", "ESSF", "ICH", "IDF", "MH", "MS", "PP", "SBPS", "SBS", "SWB")
  } else {
    zones.bc <- c("BG", "BWBS", "CDF", "CWH", "ESSF", "ICH", "IDF", "MH", "MS", "PP", "SBPS", "SBS", "SWB", "BAFA", "CMA", "IMA")
  }
  
  # colour scheme
  colScheme <- if(zone){
    zone_colours <- unique(WNA_BGCs[, .(Zone, ZoneColour)])
    zone_colours$ZoneColour
  } else {
    WNA_BGCs$SubzoneColour
  }
  
  #extract a vector of bgc labels for each cell
  bgc <- dat[,2]
  
  if(mask_alpine) bgc[grep("IMA|CMA|BAFA", bgc)] <- NA
  
  #convert to zone and factorize
  if(zone){
    bgc <- sub("^([A-Z]+).*", "\\1", bgc)
    bgc <- factor(bgc, levels = zone_colours$Zone)
  } else {
    bgc <- factor(bgc, levels = WNA_BGCs$BGC)
  }
  
  values(X) <- NA
  X[dat[,1]] <- bgc
  if(!is.null(mask)) X <- terra::mask(X, mask)
  X[1:length(levels(bgc))] <- 1:length(levels(bgc))
  
  ColScheme <- if(zone) zone_colours$ZoneColour else WNA_BGCs$SubzoneColour #(KIRI) I think this is redundant?
  
  image(X, axes=F, col=ColScheme, main = title , adj = 0.05, cex.main = 0.85, font.main = 1)

  X.mask <- X
  values(X.mask) <- NA
  values(X.mask)[1:length(levels(bgc))] <- 1 # cover up the color bar
  terra::plot(X.mask, add=T, col="white", legend=FALSE) # cover up the color bar
  
  if(!is.null(boundary)) terra::plot(boundary, add=T, border=1, lwd=0.4)
  
  if(legend){
    if(mask_alpine){
      legend("bottomleft", legend=c(zones.bc, "Alpine"), pch = 22, pt.cex = 1.5, cex = 0.9, pt.bg=c(zone_colours[match(zones.bc, zone_colours$Zone), ZoneColour], "white"), bty="n", y.intersp = 1.75, x.intersp = 1.25)
    } else {
      legend("bottomleft", legend=zones.bc, pch = 22, pt.cex = 1.5, cex = 0.9, pt.bg=zone_colours[match(zones.bc, zone_colours$Zone), ZoneColour], bty="n", y.intersp = 1.75, x.intersp = 1.25)
    }
  }
  
  if(!is.null(label_exotic)){
    exotic <- table(bgc[-which(sub("^([A-Z]+).*", "\\1", bgc)%in%zones.bc)])
    exotic <- exotic[exotic>150]
    exotic <- exotic[rev(order(exotic))]
    
    exotic.pct <- round(as.numeric(formatC(signif(exotic/length(bgc)*100,digits=3), digits=3,format="fg", flag="#")),2)
    
    bgcs.exotic <- names(exotic)
    for(bgc.exotic in bgcs.exotic){
      pts <- which(levels(bgc)[values(X)]==bgc.exotic)
      q <- q_exotic
      pt <- xyFromCell(X, pts[min(which(pts >= quantile(pts, q)))])
      points(pt, pch=21, bg=as.character(colScheme[which(levels(bgc)==bgc.exotic)]), cex=1, lwd=0.8)
      text(pt-c(0, 0), bgc.exotic, pos=4, cex=0.7, font=2, offset=0.3)
      # print(q)
    }
  }
}

#' Calculate novelty for all gcm/ssp/run/period combinations from bgc_raw
#' @description
#' Create a table in the duckdb, novelty_raw, which contains novelty values for each cell
#' @param con duckdb connection
#' @param target_pts data.table of points to calculate novelty for (should be same table used to calculate bgc_raw)
#' @param analog_pts data.table of WNA analog points 
#' @param ssps Character vector of ssps to calculate novelty for. Default is "ssp245"
#' @param append Logical. Append to existing table? Default FALSE throws error if table already exists.
#' @param table_name Character. Name of table to pull BGC predictions from. Default "bgc_raw".
#' @import data.table
#' @import climr
#' @importFrom glue glue glue_sql
#' @importFrom DBI dbGetQuery dbWriteTable
#' @export
cciss_novelty <- function(con, target_pts, analog_pts, observed = FALSE, ssps = "ssp245", append = FALSE, table_name = "bgc_raw") {
  table_name <- DBI::SQL(table_name)
  if(duckdb_table_exists(con, "novelty_raw") & !append) stop("Table novelty_raw already exists. Please drop table or set append = TRUE")
  nov_vars <- as.vector(outer(c("Tmin", "Tmax", "PPT"), c("wt", "sp", "sm", "at"), paste, sep = "_"))
  clim.pts <- downscale(xyz = analog_pts, which_refmap = "refmap_climr", return_refperiod = TRUE,
                        vars = nov_vars)
  clim.pts[analog_pts, BGC := i.BGC, on = "id"]
  
  # historical interannual climatic variability at the geographic centroids of the training points
  pts.mean <- analog_pts[, lapply(.SD, mean), by = BGC]
  pts.mean$id <- 1:dim(pts.mean)[1]
  clim.icv.pts <- downscale(xyz = pts.mean,
                            which_refmap = "refmap_climr",
                            obs_years = 1961:1990,
                            obs_ts_dataset = "cru.gpcc",
                            return_refperiod = FALSE,
                            vars = nov_vars,
                            db_option = "local")
  gcms_use <- dbGetQuery(con, glue_sql("select distinct gcm from {table_name}", .con = con))[,1]
  #ssps_use <- dbGetQuery(con, "select distinct ssp from bgc_raw")[,1]

  if(observed){
    splits <- c(seq(1, nrow(target_pts), by = 500000), nrow(target_pts) + 1)
    tmp_ls <- list()
    for(i in 1:(length(splits)-1)){
      message(glue("Downscaling chunk {i} of {length(splits)-1}"))
      tmp_ls[[i]] <- downscale(target_pts[splits[i]:(splits[i+1]-1)], 
                                which_refmap = "refmap_climr", 
                                obs_periods = "2001_2020",
                                return_refperiod = FALSE,
                                vars = nov_vars)
    }
    res <- rbindlist(tmp_ls)
    res <- res[!is.na(Tmin_sm),]
    res[is.na(res)] <- 0
    
    vars_temp <- c("id","PERIOD",nov_vars)
    ##novelty
    clim_nov <- res[,..vars_temp]
    bgc <- dbGetQuery(con, glue_sql("select * from {table_name} where period = '2001_2020_obs'", .con = con))
    setDT(bgc)
    setnames(clim_nov, old = c("id","PERIOD"), new = c("cellnum","period"))
    clim_nov[bgc, bgc_pred := i.bgc_pred, on = c("cellnum")]
    clim_nov <- na.omit(clim_nov)
    message("Calculating observed novelty...")
    clim_nov[,novelty := analog_novelty_core(clim.targets = .SD, 
                                                      clim.analogs = clim.pts, 
                                                      label.targets = bgc_pred, 
                                                      label.analogs = clim.pts$BGC, 
                                                      vars = as.vector(outer(c("Tmin", "Tmax", "PPT"), c("wt", "sp", "sm", "at"), paste, sep = "_")),
                                                      clim.icvs = clim.icv.pts,
                                                      label.icvs = pts.mean$BGC[clim.icv.pts$id],
                                                      weight.icv = 0.5,
                                                      threshold = 0.95,
                                                      pcs = NULL), by = .(period)]

    novelty_obs <- clim_nov[,.(cellnum, gcm = "Observed", ssp = "Observed", run = "Observed", period = "2001_2020_obs", novelty)]
    print(head(novelty_obs))
    dbWriteTable(con, "novelty_raw", novelty_obs, row.names = F, append = T)
  } else {
    for(gcm in gcms_use){
    for(ssp in ssps){
      periods_use <- dbGetQuery(con, glue_sql("select distinct period from {table_name} where ssp = {ssp} and gcm = {gcm}", .con = con))[,1]
      periods_use <- periods_use[!grepl("2001_2020_obs", periods_use)] # remove historical period if it exists in table
      for(period in periods_use){
        message(glue("Processing novelty for gcm = {gcm} and ssp = {ssp} and period = {period}"))
        runs <- dbGetQuery(con, glue_sql("select distinct run from {table_name} where ssp = {ssp} and gcm = {gcm} and period = {period}", .con = con))[,1]

        ## Need to save memory when downscaling - split into chunks and recombine
        splits <- c(seq(1, nrow(target_pts), by = 500000), nrow(target_pts) + 1)
        tmp_ls <- list()
        for(i in 1:(length(splits)-1)){
          message(glue("Downscaling chunk {i} of {length(splits)-1}"))
          tmp_ls[[i]] <- downscale(target_pts[splits[i]:(splits[i+1]-1)], 
                                   which_refmap = "refmap_climr", 
                                   gcms = gcm, 
                                   ssps = ssp, 
                                   gcm_periods = period, 
                                   run_nm = runs,
                                   return_refperiod = FALSE,
                                   vars = nov_vars)
        }
        res <- rbindlist(tmp_ls)
        res <- res[!is.na(Tmin_sm),]
        res[is.na(res)] <- 0
        
        vars_temp <- c("id","GCM","SSP","RUN","PERIOD",nov_vars)
        ##novelty
        clim_nov <- res[,..vars_temp]
        bgc <- dbGetQuery(con, glue_sql("select * from {table_name} where ssp = {ssp} and gcm = {gcm} and period = {period}", .con = con))
        setDT(bgc)
        setnames(clim_nov, old = c("id","GCM","SSP","RUN","PERIOD"), new = c("cellnum","gcm","ssp","run","period"))
        clim_nov[bgc, bgc_pred := i.bgc_pred, on = c("cellnum","gcm","ssp","run","period")]
        clim_nov <- na.omit(clim_nov)
        message("Calculating novelty...")
        clim_nov[,novelty := analog_novelty_core(clim.targets = .SD, 
                                                      clim.analogs = clim.pts, 
                                                      label.targets = bgc_pred, 
                                                      label.analogs = clim.pts$BGC, 
                                                      vars = as.vector(outer(c("Tmin", "Tmax", "PPT"), c("wt", "sp", "sm", "at"), paste, sep = "_")),
                                                      clim.icvs = clim.icv.pts,
                                                      label.icvs = pts.mean$BGC[clim.icv.pts$id],
                                                      weight.icv = 0.5,
                                                      threshold = 0.95,
                                                      pcs = NULL), by = .(gcm,ssp,run,period)]
        clim_nov <- clim_nov[,.(cellnum, gcm, ssp, run, period, novelty)]
        dbWriteTable(con, "novelty_raw", clim_nov, row.names = F, append = T)
      }
    }
  }
  }
  message("Written table novelty_raw to database!")
}



#' Calculate mean novelty across runs
#' @param con duckdb database connection
#' @param table_name name of table to create with mean novelty values
#' @importFrom glue glue
#' @importFrom DBI dbExecute
#' @export
ensemble_novelty <- function(con, table_name = "ensemble_novelty") {
  if(duckdb_table_exists(con, table_name)) stop(glue("Table {table_name} already exists in database. Please drop table and rerun."))
  qry <- glue("create table {table_name} as
              SELECT cellnum, period, AVG(novelty) as novelty
              FROM novelty_raw
              GROUP BY cellnum, period")
  dbExecute(con, qry)
}

#' Plot novelty map with values from database
#' @param con duckdb database connection
#' @param raster_template SpatRaster template to use for plotting
#' @param period Character. Which period to plot for?
#' @param ensemble Logical. Plot ensemble novelty (default) or raw novelty for a specified run?
#' @import terra
#' @importFrom glue glue
#' @importFrom DBI dbGetQuery
#' @importFrom data.table data.table
#' @importFrom grDevices  colorRampPalette
#' @export
plot_novelty <- function(con, raster_template, period, ensemble = TRUE, gcm = NULL, 
ssp = NULL, run = NULL, raw = FALSE, save_location = "novelty_rasters") {
  breakseq <- c(0,4,8)
  bs2 <- breakseq * 100
  breakpoints <- seq.int(bs2[1], bs2[3], 1)
  colScheme <- colorRampPalette(c("gray90", "gray50", "#FFF200", "#CD0000", "#000000"))(length(breakpoints))
  coltab <- data.table(
    values = as.integer(breakpoints),
    color = colScheme               # Corresponding colors
  )
  rt <- copy(raster_template)
  values(rt) <- NA
  if(ensemble){
    nov_curr <- dbGetQuery(con, glue("select cellnum, novelty
                                     from ensemble_novelty where period = '{period}'"))
  } else {
    nov_curr <- dbGetQuery(con, glue("select cellnum, novelty
                                     from novelty_raw where period = '{period}' and gcm = '{gcm}' and ssp = '{ssp}' and run = '{run}'"))
  }
  rt[nov_curr$cellnum] <- nov_curr$novelty
  rt <- round(rt, digits = 2)
  rt <- as.int(rt*100)
  rt[rt > 800] <- 800
  if(raw) {
    if(!dir.exists(save_location)) dir.create(save_location)
    writeRaster(rt, datatype = "INT2U", filename = paste0(save_location,"/noveltyRaw_",ifelse(ensemble, "ensemble", paste(gcm, ssp, run, sep = "_")), "_",period, ".tif"), overwrite = TRUE)
  } else {
    rgbnov <- subst(rt, coltab$values, t(col2rgb(coltab$color,alpha = TRUE)),names = c("red","green", "blue","alpha"))
    #plotRGB(rgbnov)
    if(!dir.exists(save_location)) dir.create(save_location)
    writeRaster(rgbnov, paste0(save_location,"/novelty_",ifelse(ensemble, "ensemble", paste(gcm, ssp, run, sep = "_")), "_",period, ".tif"), overwrite = TRUE)
  }
  
}