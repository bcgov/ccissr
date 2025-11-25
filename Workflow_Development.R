### CCISS Workflow Development

library(terra)
library(climr)
library(data.table)
library(ranger)
library(ccissr)
library(car)

BGCmodel <- readRDS("../Common_Files/BGCmodel_WNA_V2.1.rds")
dem <- unwrap(climr::dem_vancouver)
dem <- aggregate(dem, fact = 2)

gcms_cciss <- c("ACCESS-ESM1-5", "CNRM-ESM2-1", "EC-Earth3", "GFDL-ESM4",
                "GISS-E2-1-G", "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0")

vars_needed <- c("CMD_sm", "DDsub0_sp", "DD5_sp", "Eref_sm", "Eref_sp", "EXT", 
                 "MWMT", "NFFD_sm", "NFFD_sp", "PAS_an", "PAS_sp", "SHM", "Tave_sm", 
                 "Tave_sp", "Tmax_sm", "Tmax_sp", "Tmin_an", "Tmin_at", "Tmin_sm", 
                 "Tmin_sp", "Tmin_wt","CMI_an", "PPT_05","PPT_06","PPT_07","PPT_08",
                 "PPT_09","PPT_at","PPT_wt","CMD_07","CMD_an"
)

periods = list_gcm_periods()[-1]

summary_preds_gcm(dem, BGCmodel, vars_needed, gcms_cciss, periods_use = periods,summarise = FALSE, base_folder = base_path)

temp_ls <- list.files("cciss_spatial_vignette/bgc_data/", pattern = "bgc_raw.*", full.names = TRUE)
bgc_all_ls <- lapply(temp_ls, FUN = fread)
bgc_all <- rbindlist(bgc_all_ls)

bgc_mapped <- as.data.frame(bgc_raster$bgc_rast, cells = TRUE) 
setDT(bgc_mapped)
bgc_mapped[bgc_raster$ids, bgc := i.bgc, on = "bgc_id"]
bgc_mapped[, bgc_id := NULL]
setnames(bgc_mapped, old = "cell", new = "cellnum")
bgc_all[bgc_mapped, bgc := i.bgc, on = "cellnum"]

bgc_all <- na.omit(bgc_all)

### Spp suit
edas <- c("B2","C4","D6")
eda_table <- copy(E1)
eda_table <- eda_table[is.na(SpecialCode),]
eda_table <- eda_table[Edatopic %in% edas,]
eda_table <- unique(eda_table[,.(BGC,SS_NoSpace,Edatopic)])
setkey(eda_table, BGC)
setkey(bgc_all, bgc_pred)


bgc_eda <- merge(bgc_all, eda_table, by.x = "bgc_pred", by.y = "BGC", allow.cartesian = TRUE, all.x = TRUE)
setnames(bgc_eda, old = "SS_NoSpace", new = "SS_Pred")
bgc_eda <- merge(bgc_eda, eda_table, by.x = c("bgc", "Edatopic"), by.y = c("BGC","Edatopic"), 
                 allow.cartesian = TRUE, all.x = TRUE)
setorder(bgc_eda, Edatopic, cellnum, gcm, ssp, run, period)

suit <- copy(S1)
suit_sub <- suit[spp == "Fd",]
setkey(suit_sub,"ss_nospace")
bgc_eda[suit_sub, NewSuit := i.newfeas, on = c(SS_Pred = "ss_nospace")]
bgc_eda[suit_sub, HistSuit := i.newfeas, on = c(SS_NoSpace = "ss_nospace")]

cciss_res <- bgc_eda[,.(cellnum,ssp,gcm,run,period,Edatopic,NewSuit,HistSuit)]

if(fractional){
  cciss_res[is.na(NewSuit) | NewSuit == 4L, NewSuit := 5L]
  cciss_res[, NewSuit := 1 - (NewSuit - 1) / 4]
  
  cciss_res[is.na(HistSuit) | HistSuit == 4L, HistSuit := 5L]
  cciss_res[, HistSuit := 1 - (HistSuit - 1) / 4]
} else {
  ###convert to presence/absence
}

##persist/expand
cciss_res[,Persist := 0][HistSuit > 0, Persist := NewSuit]
cciss_res[,Expand := 0][HistSuit == 0, Expand := NewSuit]

suit_perexp <- cciss_res[,.(Persist_Tot = sum(Persist), 
                         Expand_Tot = sum(Expand)), 
                      by = .(Edatopic, ssp, gcm, run, period)]

##use original BGC map and join with species

bgc_sum <- bgc_mapped[,.(BGC_Tot = .N), by = bgc]
bgc_perexp[bgc_sum, BGC_Tot := i.BGC_Tot, on = c(bgc_pred = "bgc")]
bgc_perexp[,`:=`(Persistance = (Persist_Tot/BGC_Tot), Expansion = (Expand_Tot/BGC_Tot))]
bgc_perexp <- na.omit(bgc_perexp)

#bgc_eda[eda_table, SS_Mapped := i.SS_NoSpace, on = c(bgc = "BGC", "Edatopic")]

### Create Expansions/Persistance


bgc_all[,Persist := 0][bgc_pred == bgc, Persist := 1]
bgc_all[,Expand := 0][bgc_pred != bgc, Expand := 1]

bgc_perexp <- bgc_all[,.(Persist_Tot = sum(Persist), 
                         Expand_Tot = sum(Expand)), 
                      by = .(ssp, gcm, run, period, bgc_pred)]
bgc_sum <- bgc_mapped[,.(BGC_Tot = .N), by = bgc]
bgc_perexp[bgc_sum, BGC_Tot := i.BGC_Tot, on = c(bgc_pred = "bgc")]
bgc_perexp[,`:=`(Persistance = (Persist_Tot/BGC_Tot), Expansion = (Expand_Tot/BGC_Tot))]
bgc_perexp <- na.omit(bgc_perexp)

### Join Spp Suitability

## Can probably join expansion/persistance with spp directly? Actually probably not
## Probably need to loop through species, otherwise will be too big?

#-------------------------
# bgc bubbleplot
#-------------------------
ColScheme <- fread("../CCISS_ShinyApp/app/WNA_SZ_Cols_v13_6.csv")
unit.persistence.focal <- "none"
units <- unique(bgc_sum$bgc)
scenario <- "ssp245"
period_sel <- "2041_2060"

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

unit=units[2]
for(unit in units){
  col.focal <- if(unit.persistence.focal=="none") ColScheme$colour[which(ColScheme$classification==unit)] else "lightgray"
  col.focal2 <- if(unit.persistence.focal=="none") "black" else "darkgray"
  x <- bgc_perexp[ssp == scenario & period == period_sel & bgc_pred == unit, Persistance]
  y <- bgc_perexp[ssp == scenario & period == period_sel & bgc_pred == unit, Expansion]
  y[y<2^(ylim[1])] <- 2^(ylim[1])
  y <- log2(y)
  
  # points(x,y)
  if(length(x)>1 & var(x) > 0){
    if(var(y)==0) {
      lines(range(x, na.rm=T), range(y), col=col.focal)
    } else {
      print(paste0("Data Ellipse for ", unit))
      dataEllipse(x, y, levels=0.5, center.pch=21, add=T, col=col.focal, fill=T, lwd=0.5, plot.points=F)
    } 
  }
  points(mean(x),mean(y), pch=21, bg=col.focal, cex=if(unit==unit.persistence.focal) 4.5 else 3, col=col.focal2)
  text(mean(x),mean(y), unit, cex=if(unit==unit.persistence.focal) 1 else 0.7, font=2, col=col.focal2)
}

if(unit.persistence.focal!="none"){
  
  for(gcm in gcms){
    i=which(gcms==gcm)
    runs <- unique(identity$RUN[which(identity$GCM==gcm)])
    for(run in runs){
      is.focalSim <- paste(gcm, substr(run,1,2), sep="_")==sim.focal
      x2 <- persistence[c(1, which(identity$GCM==gcm & identity$SSP==scenario & identity$RUN==run)), which(names(persistence)==unit.persistence.focal)]
      y2 <- expansion[c(1, which(identity$GCM==gcm & identity$SSP==scenario & identity$RUN==run)), which(names(expansion)==unit.persistence.focal)]
      y2[y2<2^(ylim[1])] <- 2^(ylim[1])
      y2 <- log2(y2)
      # if(all(diff(x2) > 0)){
      if(length(unique(sign(diff(x2))))==1 & sum(diff(x2))!=0){
        x3 <- if(unique(sign(diff(x2)))==-1) rev(x2) else x2
        y3 <- if(unique(sign(diff(x2)))==-1) rev(y2) else y2
        s <- stinterp(x3,y3, seq(min(x3),max(x3), diff(xlim)/500)) # way better than interpSpline, not prone to oscillations
        if(run=="ensembleMean") lines(s, col=ColScheme.gcms[i], lwd=2)
      } else if(run=="ensembleMean") lines(x2, y2, col=ColScheme.gcms[i], lwd=2)
      j=which(periods==period)
      if(length(runs)>2) points(x2[j],y2[j], pch=21, bg=ColScheme.gcms[i], cex=1)
      if(run=="ensembleMean") points(x2[j],y2[j], pch=21, bg=ColScheme.gcms[i], cex=if(gcm=="Ensemble mean") 3.5 else 3)
      if(run=="ensembleMean") text(x2[j],y2[j], mods[i], cex=if(gcm=="Ensemble mean") 0.7 else 0.5, font=2)
    }
  }
  
  if(recent==T){
    x1 <- persistence[2, which(names(persistence)==unit.persistence.focal)]
    y1 <- expansion[2, which(names(expansion)==unit.persistence.focal)]
    y1[y1<2^(ylim[1])] <- 2^(ylim[1])
    y1 <- log2(y1)
    points(x1,y1, pch=21, bg="gray", col=1, cex=2)
    text(x1,y1, "2001-2020 (observed)", cex=1, font=2, pos=4, col="gray30", offset=0.7)
  }
  
}
box()





## Colin Spp Bubbleplots
## ---------------------

if(edatope=="B2") spp.focal <- input$spp.focal.2.B2
if(edatope=="C4") spp.focal <- input$spp.focal.2.C4
if(edatope=="D6") spp.focal <- input$spp.focal.2.D6

persistence <- if(fractional==T) get(paste("suit.persistence", edatope, sep=".")) else get(paste("spp.persistence", edatope, sep="."))
expansion <- if(fractional==T) get(paste("suit.expansion", edatope, sep=".")) else get(paste("spp.expansion", edatope, sep="."))
spps <- names(persistence)

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
# arctext(x = "Growing feasible range", center = c(-1, -28.7), radius = 4.6, start = 0.431*pi , cex = 0.8, stretch = 1.05, col="darkgray", font=2)
# arctext(x = "Shrinking feasible range", center = c(-1, -29.3), radius = 4.6, start = 0.431*pi , cex = 0.8, stretch = 1.05, col="darkgray", font=2)
# mtext(paste(edatope.names[which(edatopes==edatope)], " sites", " (", edatope, ")", sep=""), side=3, line=-1.25, adj= if(edatope=="C4") 0.025 else 0.075, cex=0.7, font=1)

spp=spps[1]
for(spp in spps){
  i <- which(spps==spp)
  col.focal <- if(spp.focal=="none") sppcolors[i] else "lightgray"
  col.focal2 <- if(spp.focal=="none") "black" else "darkgray"
  x <- persistence[which(identity$SSP==scenario & identity$PERIOD==period), i]
  y <- expansion[which(identity$SSP==scenario & identity$PERIOD==period), i]
  y[y<2^(ylim[1])] <- 2^(ylim[1])
  y <- log2(y)
  
  # points(x,y)
  if(length(x)>1){
    if(var(x)>0) if(var(y)==0) lines(range(x), range(y), col=col.focal) else dataEllipse(x, y, levels=0.5, center.pch=21, add=T, col=col.focal, fill=T, lwd=0.5, plot.points=F)
  }
  points(mean(x),mean(y), pch=21, bg=col.focal, cex=if(spp==spp.focal) 4.5 else 3, col=col.focal2)
  text(mean(x),mean(y), spp, cex=if(spp==spp.focal) 1 else 0.7, font=2, col=col.focal2)
}

if(spp.focal!="none"){
  for(gcm in gcms){
    i=which(gcms==gcm)
    runs <- unique(identity$RUN[which(identity$GCM==gcm)])
    for(run in runs){
      is.focalSim <- paste(gcm, substr(run,1,2), sep="_")==sim.focal
      x2 <- persistence[c(1, which(identity$GCM==gcm & identity$SSP==scenario & identity$RUN==run)), which(spps==spp.focal)]
      y2 <- expansion[c(1, which(identity$GCM==gcm & identity$SSP==scenario & identity$RUN==run)), which(spps==spp.focal)]
      y2[y2<2^(ylim[1])] <- 2^(ylim[1])
      y2 <- log2(y2)
      # if(all(diff(x2) > 0)){
      if(length(unique(sign(diff(x2))))==1 & sum(diff(x2))!=0){
        x3 <- if(unique(sign(diff(x2)))==-1) rev(x2) else x2
        y3 <- if(unique(sign(diff(x2)))==-1) rev(y2) else y2
        s <- stinterp(x3,y3, seq(min(x3),max(x3), diff(xlim)/500)) # way better than interpSpline, not prone to oscillations
        if(run=="ensembleMean") lines(s, col=ColScheme.gcms[i], lwd=2)
      } else if(run=="ensembleMean") lines(x2, y2, col=ColScheme.gcms[i], lwd=2)
      j=which(periods==period)
      if(length(runs)>2) points(x2[j],y2[j], pch=21, bg=ColScheme.gcms[i], cex=1)
      if(run=="ensembleMean") points(x2[j],y2[j], pch=21, bg=ColScheme.gcms[i], cex=if(gcm=="ensembleMean") 3.5 else 3)
      if(run=="ensembleMean") text(x2[j],y2[j], mods[i], cex=if(gcm=="ensembleMean") 0.7 else 0.5, font=2)
    }
  }
}

