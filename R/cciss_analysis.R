### Kiri Daust, 2025

#' Create table of species persistance/expansion from raw BGC projections
#' @param edatopes Character. Vector of desired edatopic positions.
#' @param fractional Logical. Use fractional (suitability based) values for calculations?
#' @param bgc_mapped List containing SpatRaster of BGCs and id table, or data.table (must have columns `cell`,`BGC`). Usually created using `make_bgc_template`
#' @param species Character. Either "auto", or a vector of species codes to use. If "auto", all species with historically suitable range will be used.
#' @param periods Character. Either "auto" or vector of period names. If "auto", all periods present in bgc projections will be used.
#' @param base_folder Base folder to read results from.
#' @return data.table containing persistance and expansion for each model/period/scenario/species
#' @import data.table
#' @export
spp_persist_expand <- function(edatopes, fractional, bgc_mapped, species = "auto", periods = "auto", base_folder){
  in_folder <- file.path(base_folder,"bgc_data")
  obs <- FALSE
  ##read in BGC preds
  if(obs){
    periods <- list_obs_periods()
    bgc_all <- fread(file.path(in_folder, "bgc_summary_obs.csv"))
    obs_nm <- "obs_"
  } else {
    temp_ls <- list.files(in_folder, pattern = "bgc_raw.*", full.names = TRUE)
    #temp_ls <- temp_ls[!grepl("obs",temp_ls)]
    bgc_all_ls <- lapply(temp_ls, FUN = fread)
    bgc_all <- rbindlist(bgc_all_ls)
    rm(bgc_all_ls)
    obs_nm <- ""
  }
  
  if(periods[1] != "auto") {
    bgc_all <- bgc_all[period %in% periods,]
  }
  
  ##create table of mapped BGCs
  if(inherits(bgc_mapped, "list")){
    bgc_rast <- bgc_mapped$bgc_rast
    rast_ids <- bgc_mapped$ids
    bgc_points <- as.data.frame(bgc_rast, cells=T) |> as.data.table()
    bgc_points[rast_ids, bgc := i.bgc, on = "bgc_id"]
    bgc_points[,bgc_id := NULL]
  } else {
    bgc_points <- bgc_mapped |> as.data.table()
  }
  
  setnames(bgc_points, old = "cell", new = "cellnum")
  bgc_all[bgc_points, bgc := i.bgc, on = "cellnum"]
  
  bgc_all <- na.omit(bgc_all)
  
  ### Prep edatopic and suitability tables
  eda_table <- copy(E1) ##Edatopic table
  eda_table <- eda_table[is.na(SpecialCode),]
  eda_table <- eda_table[Edatopic %in% edatopes,]
  eda_table <- unique(eda_table[,.(BGC,SS_NoSpace,Edatopic)])
  setkey(eda_table, BGC)
  setkey(bgc_all, bgc_pred)
  
  suit <- copy(S1) ##Suitability table
  #suit <- suit[spp %in% species,]
  suit <- na.omit(suit, cols = "spp")
  
  #########################################
  bgc_eda <- merge(bgc_all, eda_table, by.x = "bgc_pred", by.y = "BGC", allow.cartesian = TRUE, all.x = TRUE)
  setnames(bgc_eda, old = "SS_NoSpace", new = "SS_Pred")
  bgc_eda <- merge(bgc_eda, eda_table, by.x = c("bgc", "Edatopic"), by.y = c("BGC","Edatopic"), 
                   allow.cartesian = TRUE, all.x = TRUE)
  setorder(bgc_eda, Edatopic, cellnum, gcm, ssp, run, period)
  
  ##mapped suitability
  bgc_sum <- bgc_points[,.(BGC_Tot = .N), by = bgc]
  mapped_ss <- merge(bgc_sum, eda_table, by.x = "bgc", by.y = "BGC", allow.cartesian = TRUE)
  mapped_ss <- merge(mapped_ss, suit[,.(ss_nospace,spp,newfeas)], by.x = "SS_NoSpace", by.y = "ss_nospace", 
                     all.x = TRUE, allow.cartesian = TRUE)
  setnames(mapped_ss, old = "newfeas", new = "MappedSuit")
  mapped_ss[is.na(MappedSuit) | MappedSuit == 4, MappedSuit := 5]
  mapped_suit <- mapped_ss[,.(Suit = min(MappedSuit)), by = .(spp, bgc, BGC_Tot, Edatopic)]
  
  if(species[1] == "auto"){
    spp_list <- unique(mapped_suit[Suit <= 3, spp])
  } else {
    spp_list <- species
  }
  
  if(fractional){
    mapped_suit[,Suit := 1 - (Suit - 1) / 4]
    mapped_suit[,FracSuit := Suit * BGC_Tot]
    mapped_suit <- mapped_suit[,.(MappedSuit = sum(FracSuit)), by = .(spp, Edatopic)]
  } else {
    mapped_suit[Suit != 5,Suit := 1]
    mapped_suit[Suit == 5,Suit := 0]
    mapped_suit[,FracSuit := Suit * BGC_Tot]
    mapped_suit <- mapped_suit[,.(MappedSuit = sum(FracSuit)), by = .(spp, Edatopic)]
  }
  
  suit <- suit[spp %in% spp_list,] ##subset for efficiency
  
  spp_res <- list()
  for(spp_curr in spp_list){
    message(spp_curr)
    suit_sub <- suit[spp == spp_curr,]
    setkey(suit_sub,"ss_nospace")
    if("NewSuit" %in% names(bgc_eda)){
      bgc_eda[,c("NewSuit","HistSuit") := NULL]
    }
    bgc_eda[suit_sub, NewSuit := i.newfeas, on = c(SS_Pred = "ss_nospace")]
    bgc_eda[suit_sub, HistSuit := i.newfeas, on = c(SS_NoSpace = "ss_nospace")]
    
    cciss_res <- bgc_eda[,.(cellnum,ssp,gcm,run,period,Edatopic,NewSuit,HistSuit)]
    
    if(fractional){
      cciss_res[is.na(NewSuit) | NewSuit == 4L, NewSuit := 5L]
      cciss_res[, NewSuit := 1 - (NewSuit - 1) / 4]
      
      cciss_res[is.na(HistSuit) | HistSuit == 4L, HistSuit := 5L]
      cciss_res[, HistSuit := 1 - (HistSuit - 1) / 4]
    } else {
      cciss_res[is.na(NewSuit) | NewSuit == 4L, NewSuit := 5L][
        NewSuit != 5, NewSuit := 1
      ][
        NewSuit == 5, NewSuit := 0
      ]
      cciss_res[is.na(HistSuit) | HistSuit == 4L, HistSuit := 5L][
        HistSuit != 5, HistSuit := 1
      ][
        HistSuit == 5, HistSuit := 0
      ]
    }
    
    ## best suitability for each location (often multiple overlapping site series)
    cciss_res <- cciss_res[,.(NewSuit = max(NewSuit), HistSuit = max(HistSuit)), 
                           by = .(cellnum, ssp, gcm, run, period, Edatopic)]
    ##persist/expand
    cciss_res[,Persist := 0][HistSuit > 0, Persist := NewSuit]
    cciss_res[,Expand := 0][HistSuit == 0, Expand := NewSuit]
    
    suit_perexp <- cciss_res[,.(Persist_Tot = sum(Persist), 
                                Expand_Tot = sum(Expand)), 
                             by = .(Edatopic, ssp, gcm, run, period)]
    
    suit_perexp[mapped_suit[spp == spp_curr,], Mapped := i.MappedSuit, on = "Edatopic"]
    suit_perexp[,`:=`(Persistance = Persist_Tot/Mapped, Expansion = Expand_Tot/Mapped)]
    suit_perexp[,c("Persist_Tot","Expand_Tot","Mapped") := NULL]
    
    spp_res[[spp_curr]] <- suit_perexp
  }
  
  spp_perexp <- rbindlist(spp_res, idcol = "spp")
  spp_perexp <- na.omit(spp_perexp)
  return(spp_perexp)
}

#' Create table of BGC subzone/variant persistance/expansion from raw BGC projections
#' @param bgc_mapped List containing SpatRaster of BGCs and id table, or data.table (must have columns `cell`,`BGC`). Usually created using `make_bgc_template`
#' @param periods Character. Either "auto" or vector of period names. If "auto", all periods present in bgc projections will be used.
#' @param base_folder Base folder to read results from.
#' @return data.table containing persistance and expansion for each model/period/scenario/bgc
#' @import data.table
#' @export
bgc_persist_expand <- function(bgc_mapped, periods = "auto", base_folder){
  in_folder <- file.path(base_folder,"bgc_data")
  obs <- FALSE
  ##read in BGC preds
  if(obs){
    periods <- list_obs_periods()
    bgc_all <- fread(file.path(in_folder, "bgc_summary_obs.csv"))
    obs_nm <- "obs_"
  } else {
    temp_ls <- list.files(in_folder, pattern = "bgc_raw.*", full.names = TRUE)
    #temp_ls <- temp_ls[!grepl("obs",temp_ls)]
    bgc_all_ls <- lapply(temp_ls, FUN = fread)
    bgc_all <- rbindlist(bgc_all_ls)
    rm(bgc_all_ls)
    obs_nm <- ""
  }
  
  if(periods[1] != "auto") {
    bgc_all <- bgc_all[period %in% periods,]
  }
  
  ##create table of mapped BGCs
  if(inherits(bgc_mapped, "list")){
    bgc_rast <- bgc_mapped$bgc_rast
    rast_ids <- bgc_mapped$ids
    bgc_points <- as.data.frame(bgc_rast, cells=T) |> as.data.table()
    bgc_points[rast_ids, bgc := i.bgc, on = "bgc_id"]
    bgc_points[,bgc_id := NULL]
  } else {
    bgc_points <- bgc_mapped |> as.data.table()
  }
  
  setnames(bgc_points, old = "cell", new = "cellnum")
  bgc_all[bgc_points, bgc := i.bgc, on = "cellnum"]
  
  bgc_all <- na.omit(bgc_all)
  
  bgc_all[,Persist := 0][bgc_pred == bgc, Persist := 1]
  bgc_all[,Expand := 0][bgc_pred != bgc, Expand := 1]
  
  bgc_perexp <- bgc_all[,.(Persist_Tot = sum(Persist), 
                           Expand_Tot = sum(Expand)), 
                        by = .(ssp, gcm, run, period, bgc_pred)]
  bgc_sum <- bgc_points[,.(BGC_Tot = .N), by = bgc]
  bgc_perexp[bgc_sum, BGC_Tot := i.BGC_Tot, on = c(bgc_pred = "bgc")]
  bgc_perexp[,`:=`(Persistance = (Persist_Tot/BGC_Tot), Expansion = (Expand_Tot/BGC_Tot))]
  bgc_perexp <- na.omit(bgc_perexp)
  
  return(bgc_perexp)
}

#' Create table of species relative suitabile area
#' @param edatopes Character. Vector of desired edatopic positions.
#' @param fractional Logical. Use fractional (suitability based) values for calculations?
#' @param spp_list Character vector of species to analyse
#' @param bgc_mapped List containing SpatRaster of BGCs and id table, or data.table (must have columns `cell`,`BGC`). Usually created using `make_bgc_template`
#' @param base_folder Base folder to read results from.
#' @return data.table containing relative area for each run/model/period/scenario/species
#' @import data.table
#' @export
spp_suit_area <- function(edatopes, fractional, spp_list, bgc_mapped, base_folder) {
  in_folder <- file.path(base_folder,"bgc_data")
  obs <- FALSE
  ##read in BGC preds
  if(obs){
    periods <- list_obs_periods()
    bgc_all <- fread(file.path(in_folder, "bgc_summary_obs.csv"))
    obs_nm <- "obs_"
  } else {
    temp_ls <- list.files(in_folder, pattern = "bgc_raw.*", full.names = TRUE)
    #temp_ls <- temp_ls[!grepl("obs",temp_ls)]
    bgc_all_ls <- lapply(temp_ls, FUN = fread)
    bgc_all <- rbindlist(bgc_all_ls)
    rm(bgc_all_ls)
    obs_nm <- ""
    bgc_all <- bgc_all[run != "ensembleMean",]
  }
  
  ##create table of mapped BGCs
  if(inherits(bgc_mapped, "list")){
    bgc_rast <- bgc_mapped$bgc_rast
    rast_ids <- bgc_mapped$ids
    bgc_points <- as.data.frame(bgc_rast, cells=T) |> as.data.table()
    bgc_points[rast_ids, bgc := i.bgc, on = "bgc_id"]
    bgc_points[,bgc_id := NULL]
  } else {
    bgc_points <- bgc_mapped |> as.data.table()
  }
  
  setnames(bgc_points, old = "cell", new = "cellnum")
  bgc_all[bgc_points, bgc := i.bgc, on = "cellnum"]
  
  bgc_all <- na.omit(bgc_all)
  
  ### Prep edatopic and suitability tables
  eda_table <- copy(E1) ##Edatopic table
  eda_table <- eda_table[is.na(SpecialCode),]
  eda_table <- eda_table[Edatopic %in% edatopes,]
  eda_table <- unique(eda_table[,.(BGC,SS_NoSpace,Edatopic)])
  setkey(eda_table, BGC)
  setkey(bgc_all, bgc_pred)
  
  suit <- copy(S1) ##Suitability table
  suit <- suit[spp %in% spp_list,]
  suit <- na.omit(suit, cols = "spp")
  
  bgc_eda <- merge(bgc_all, eda_table, by.x = "bgc_pred", by.y = "BGC", allow.cartesian = TRUE, all.x = TRUE)
  setnames(bgc_eda, old = "SS_NoSpace", new = "SS_Pred")
  
  ##mapped suitability
  bgc_sum <- bgc_points[,.(BGC_Tot = .N), by = bgc]
  mapped_ss <- merge(bgc_sum, eda_table, by.x = "bgc", by.y = "BGC", allow.cartesian = TRUE)
  mapped_ss <- merge(mapped_ss, suit[,.(ss_nospace,spp,newfeas)], by.x = "SS_NoSpace", by.y = "ss_nospace", 
                     all.x = TRUE, allow.cartesian = TRUE)
  setnames(mapped_ss, old = "newfeas", new = "MappedSuit")
  mapped_ss[is.na(MappedSuit) | MappedSuit == 4, MappedSuit := 5]
  mapped_suit <- mapped_ss[,.(Suit = min(MappedSuit)), by = .(spp, bgc, BGC_Tot, Edatopic)]
  
  if(fractional){
    mapped_suit[,Suit := 1 - (Suit - 1) / 4]
    mapped_suit[,FracSuit := Suit * BGC_Tot]
    mapped_suit <- mapped_suit[,.(MappedSuit = sum(FracSuit)), by = .(spp, Edatopic)]
  } else {
    mapped_suit[Suit != 5,Suit := 1]
    mapped_suit[Suit == 5,Suit := 0]
    mapped_suit[,FracSuit := Suit * BGC_Tot]
    mapped_suit <- mapped_suit[,.(MappedSuit = sum(FracSuit)), by = .(spp, Edatopic)]
  }
  mapped_suit <- na.omit(mapped_suit)
  
  spp_res <- list()
  for(spp_curr in spp_list){
    message(spp_curr)
    suit_sub <- suit[spp == spp_curr,]
    setkey(suit_sub,"ss_nospace")
    if("NewSuit" %in% names(bgc_eda)) bgc_eda[, NewSuit := NULL]
    
    bgc_eda[suit_sub, NewSuit := i.newfeas, on = c(SS_Pred = "ss_nospace")]
    if(fractional){
      bgc_eda[is.na(NewSuit) | NewSuit == 4L, NewSuit := 5L]
      bgc_eda[, NewSuit := 1 - (NewSuit - 1) / 4]
    } else {
      bgc_eda[is.na(NewSuit) | NewSuit == 4L, NewSuit := 5L][
        NewSuit != 5, NewSuit := 1
      ][
        NewSuit == 5, NewSuit := 0
      ]
    }
    
    cciss_res <- bgc_eda[,.(NewSuit = max(NewSuit)), 
                         by = .(cellnum, ssp, gcm, run, period, Edatopic)]
    
    suit_area <- cciss_res[,.(Proj_Area = sum(NewSuit)), by = .(Edatopic, ssp, gcm, run, period)]
    suit_area[mapped_suit[spp == spp_curr,], Mapped := i.MappedSuit, on = "Edatopic"]
    suit_area[,Suit_Prop := Proj_Area/Mapped]
    spp_res[[spp_curr]] <- suit_area
  }
  res <- rbindlist(spp_res, idcol = "spp")
  res <- na.omit(res)
  return(res)
}

#' Create a spaghetti plot for a single species over multiple edatopes
#' @param suit_area data.table. Usually created with `spp_suit_area`
#' @param species Character. Single species to create plot for
#' @return NULL. Creates plot
#' @import data.table
#' @import ggplot2
#' @importFrom stinepack stinterp
#' @export
spp_spaghettiplot <- function(suit_area, species) {
  suit_area <- suit_area[spp == species,]
  suit_area[,Year := as.integer(substr(period,1,4))]
  dat_spline <- suit_area[
    , {
      dt <- rbind(.SD, data.table(Year = 2000, Suit_Prop = 1), fill = TRUE)
      o <- order(dt$Year)
      x <- dt$Year[o]
      y <- dt$Suit_Prop[o]
      
      # interpolate
      xout <- seq(min(x), max(x), length.out = 100)
      yout <- stinterp(x, y, xout)$y
      
      .(Year = xout, Suit_Spline = yout)
    },
    by = .(Edatopic, ssp, gcm, run)
  ]
  
  dat_spline[, Group := interaction(ssp, gcm, run)]
  mean_spline <- dat_spline[
    , .(Suit_Spline = mean(Suit_Spline)),
    by = .(Edatopic, ssp, Year)
  ]
  
  ggplot(dat_spline, aes(x = Year, y = Suit_Spline, group = Group, colour = ssp)) +
    geom_line(alpha = 0.4) +
    geom_line(data = mean_spline, aes(group = ssp), size = 1.4) +
    facet_wrap(~Edatopic)+
    theme_minimal() +
    ylab("Proportion of Historic Suitable Area")
  
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
  ColScheme <- copy(subzones_colours_ref)
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
