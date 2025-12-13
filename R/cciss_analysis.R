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
  ColScheme <- rbind(copy(subzones_colours_ref),data.table(classification = names(colScheme), colour = colScheme))
  
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
#' @import data.table
#' @import ggplot2
#' @importFrom ggalluvial geom_alluvium
#' @export
plot_spparea <- function(spp, edatope, bgc_template, fractional, by_zone = TRUE, base_folder) {
  in_folder <- file.path(base_folder,"cciss_suit")
  temp_ls <- list.files(in_folder, pattern = paste0("CCISS_.*",edatope,".csv"), full.names = TRUE)
  cciss_suit <- lapply(temp_ls, FUN = fread)
  cciss_suit <- rbindlist(cciss_suit)
  cciss_spp <- cciss_suit[Spp == spp,]
  
  bgc_mapped <- as.data.frame(bgc_template$bgc_rast, cells = TRUE) 
  setDT(bgc_mapped)
  bgc_mapped[bgc_template$ids, bgc := i.bgc, on = "bgc_id"]
  bgc_mapped[, bgc_id := NULL]
  if(by_zone) {
    bgc_mapped[, zone := regmatches(bgc, regexpr("^[A-Z]+", bgc))]
  } else {
    bgc_mapped[, zone := bgc]
  }
  
  cciss_spp[bgc_mapped, zone := i.zone, on = c(SiteRef = "cell")]
  
  
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
  
  # Suppose these are your time points:
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
  

  
  cellarea <- (res(bgc_template$bgc_rast)[2]*111)*(res(bgc_template$bgc_rast)[1]*111*cos(mean(ext(bgc_template$bgc_rast)[3:4]) * pi / 180))
  cciss_sum[, SppArea := SppArea * cellarea]
  
  ggplot(cciss_sum, aes(x = Year, y = SppArea, fill = zone)) +
    geom_alluvium(aes(alluvium = zone), alpha= .9, color = "black") +
    geom_rect(
      data = bars,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "white", alpha = 0.7    # adjust alpha for translucency
    ) +
    theme_bw() +
    scale_fill_manual(values = colScheme) +
    ylab("Species Suitable Area (Km^2)")
  
}

#' Create C.R. Mahony's two-panel map plot of historic suitability, suitability change, and a boxplot of change by zone
#' @description
#' Note that currently, this function only works correctly if the preceeding analysis has been done using an Albers grid.
#' 
#' @param spp Character. Species code to create plot for
#' @param edatope Character. Edatope code to create plot for (e.g., "C4")
#' @param period Character. Period to create plot for.
#' @param bgc_template List containing SpatRaster of BGCs and id table. Usually created using `make_bgc_template`
#' @param outline SpatVector of aoi boundary.
#' @param save_png Logical. Save plot to png? Default `TRUE`. If `FALSE`, creates plot on default plotting device.
#' @param base_folder Base folder to read results from.
#' @return NULL. Creates plot
#' @import data.table
#' @import terra
#' @export
plot_2panel <- function(spp, edatope, period, bgc_template, outline, save_png = TRUE, base_path) {
  
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
    png(file=paste("./Two_Panel",spp,edatope,period,"png",sep = "."), type="cairo", units="in", width=6.5, height=5, pointsize=12, res=300)
  }
  
  par(plt=c(0,1,0,1), bg="white")
  plot(0, col="white", xaxt="n", yaxt="n", xlab="", ylab="")
  Common <- as.character(spps.lookup$EnglishName[which(spps.lookup$TreeCode==spp)])
  Latin <- as.character(spps.lookup$ScientificName[which(spps.lookup$TreeCode==spp)])
  mtext(if(spp%in%spps.lookup$TreeCode) bquote(bold(.(spp))~"-"~.(Common)) else bquote(bold(.(spp))),
        side=3, line=-2.5, adj=0.01, cex=0.9, font=2)
  if(edatope %in% edatope.names) {
    mtext(paste("Site type: ", edatope, " (", edatope.names[edatopes == edatope], ")", sep=""), side=3, line=-3.5, adj=0.01, cex=0.8, font=1)
  } else {
    mtext(paste("Site type: ", edatope, sep=""), side=3, line=-3.5, adj=0.01, cex=0.8, font=1)
  }
  mtext(paste("Time period: ", period, sep=""), side=3, line=-4.5, adj=0.01, cex=0.8, font=1)
  
  ##=================================
  ###historic suitability
  feasVals <- fread(paste0(base_path,"/cciss_suit/CCISS_", period, "_",edatope,".csv"))
  dat_spp <- feasVals[Spp == spp,]
  dat_spp[,FeasChange := Curr - Newsuit]
  dat_spp[,Curr := as.integer(round(Curr))]
  dat_spp[is.na(Curr) | Curr > 3.5, Curr := 5]
  X <- copy(bgc_template$bgc_rast)
  values(X) <- NA

  X[dat_spp$SiteRef] <- dat_spp$Curr
  breakseq <- c(0.5,1.5,2.5,3.5,5)
  ColScheme <- c("darkgreen", "dodgerblue1", "gold2", "white")
  
  par(plt = c(0, 0.5, 0.05, 0.6),new = TRUE, xpd = TRUE)
  
  image(X,xlab = NA,ylab = NA,bty = "n",  xaxt="n", yaxt="n",
        col=ColScheme, breaks=breakseq,asp = 1)
  terra::plot(outline, add=T, border="black",col = NA, lwd=0.4)
  legend("topleft", legend=c("1 (primary)", "2 (secondary)", "3 (tertiary)"),
         fill=ColScheme, bty="n", cex=0.8, title="Historical feasibility", inset=c(0.1,-0.3))
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
  
  par(plt = c(0.25, 0.95, 0.175, 1), xpd = TRUE, new = TRUE)
  image(X,xlab = NA,ylab = NA,bty = "n", xaxt="n", yaxt="n", col=ColScheme, breaks=breakpoints, asp = 1)
  image(X2, add=T, xlab = NA,ylab = NA,bty = "n", xaxt="n", yaxt="n", col=ColScheme2, breaks=breakpoints, asp = 1)
  image(X3, add=T, xlab = NA,ylab = NA,bty = "n", xaxt="n", yaxt="n", col=ColScheme3, asp = 1)
  terra::plot(outline, add=T, border="black",col = NA, lwd=0.4)
  
  xl <- 1600000; yb <- 1000000; xr <- 1700000; yt <- 1700000; xadj <- 10000
  y.int <- (yt-yb)/length(ColScheme)
  rect(xl+xadj,  head(seq(yb,yt,y.int),-1),  xr,  tail(seq(yb,yt,y.int),-1),  col=ColScheme)
  rect(xl-diff(c(xl+xadj, xr)),  head(seq(yb,yt,y.int),-1),  xl-xadj,  tail(seq(yb,yt,y.int),-1),  col=ColScheme2)
  rect(xl-diff(c(xl+xadj, xr)),  yb,  xl-xadj,  (yb+yt)/2,  col="white")
  text(xl-diff(c(xl+xadj, xr))/2, yb+(yt-yb)/4, "Expansion", srt=90, cex=0.85, font=1)
  text(rep(xr-10000,length(labels)),seq(yb,yt,(yt-yb)/(length(labels)-1)),labels,pos=4,cex=0.8,font=1)
  text(xl-diff(c(xl+xadj, xr))-30000, mean(c(yb,yt))-30000, paste("Mean change in feasibility", sep=""), srt=90, pos=3, cex=0.85, font=2)
  rect(xl+xadj,  yb-y.int-20000,  xr,  yb-20000,  col="black")
  text(xr, yb-y.int/2-30000, "Loss", pos=4, cex=0.8, font=1)
  par(xpd=F)

  ##=================================
  ## Summary by zone
  bgc_mapped <- as.data.frame(bgc_template$bgc_rast, cells = TRUE) 
  setDT(bgc_mapped)
  bgc_mapped[bgc_template$ids, bgc := i.bgc, on = "bgc_id"]
  bgc_mapped[, bgc_id := NULL]
  bgc_mapped[, zone := regmatches(bgc, regexpr("^[A-Z]+", bgc))]
  
  dat_spp[bgc_mapped, zone := i.zone, on = c(SiteRef = "cell")]
  
  
  par(mar=c(0,0,0,0), plt = c(0.77, 0.995, 0.001, 0.31), new = TRUE, mgp=c(1.25,0.15,0))
  plot(0, xlim=c(0,1), ylim=c(0,1), col="white", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
  
  par(mar=c(4.5,2,0.1,0.1), plt = c(0.7, 0.995, 0.08, 0.2), new = TRUE, mgp=c(1.25,0.15,0))
  ylim=c(-3,3)
  zones_curr <- unique(dat_spp$zone)
  xlim=c(1, length(unique(dat_spp$zone)))
  z <- boxplot(FeasChange~zone, data = dat_spp, ylab="", vertical = TRUE, plot=F)
  for(i in 1:length(zones_curr)){
    temp <- dat_spp[zone == zones_curr[i], FeasChange]
    z$stats[c(1,5), i] <- quantile(temp[!is.na(temp)],c(0.05, 0.95))
  }
  bxp(z, xlim=xlim, ylim=ylim, xaxt="n", yaxt="n", xaxs="i", ylab="", pch=0,outline=FALSE)
  lines(c(-99,99), c(0,0), lwd=2, col="darkgrey")
  bxp(z, add=T, boxfill = zoneScheme[match(zones_curr, names(zoneScheme))], xaxt="n", yaxt="n", xaxs="i", ylab="", pch=0,outline=FALSE)
  axis(1, at=1:length(zones_curr), zones_curr, tick=F, las=2, cex.axis=0.65)
  axis(2,at=seq(ylim[1], ylim[2], 3), seq(ylim[1], ylim[2], 3), las=2, tck=0)
  mtext("Mean change in feasibility", side=3, line=0.1, adj=.975, cex=0.65, font=2)
  
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
