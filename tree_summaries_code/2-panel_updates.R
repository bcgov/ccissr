# Jan 8 2026
# Mia modifying the 2-panel plot function for the following 
# goals -----------------------
#1) Fix the error I've been getting (done)
#2) Plot the reference period maps separately (STUCK)
#3) Modifications to the future-period map + boxplot: 
  # I think it would be nice to have the E1-3 suitability map (left) shown for the same time period as the mean change map (right). I.e., keep the layout of the figure as-is but with a consistent time period. (Update: I decided against this for now)
  # ALT: remove the historical suitability map and make the boxplot larger on the left panel instead (done)
  #a.	Could we also change the axes throughout this plot to “suitability” from “feasibility”? (done)

 




# 1) Modifications to Kiri's function: ----------------------

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
  
  terra::image(X,xlab = NA,ylab = NA,bty = "n",  xaxt="n", yaxt="n",
        col=ColScheme, breaks=breakseq,asp = 1)
  terra::plot(outline, add=T, border="black",col = NA, lwd=0.4)
  par(xpd = NA)
  legend("topleft", legend=c("1 (primary)", "2 (secondary)", "3 (tertiary)"), fill=ColScheme, bty="n", cex=0.8, title="Historical suitability", inset=c(0,-0.3))
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
  
  terra::image(X,xlab = NA,ylab = NA,bty = "n", xaxt="n", yaxt="n", col=ColScheme, breaks=breakpoints, asp = 1)
  terra::image(X2, add=T, xlab = NA,ylab = NA,bty = "n", xaxt="n", yaxt="n", col=ColScheme2, breaks=breakpoints, asp = 1)
  terra::image(X3, add=T, xlab = NA,ylab = NA,bty = "n", xaxt="n", yaxt="n", col=ColScheme3, asp = 1)
  terra::plot(outline, add=T, border="black",col = NA, lwd=0.4)
  
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
  text(xl-diff(c(xl+xadj, xr))-30000, mean(c(yb,yt))-30000, paste("Mean change in suitability", sep=""), srt=90, pos=3, cex=0.85, font=2)
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
  mtext("Mean change in suitability", side=3, line=0.1, adj=.975, cex=0.65, font=2)
  
  if(three_panel) {
    values(X) <- NA
    X[dat_spp[FeasChange > 0,SiteRef]] <- dat_spp[FeasChange > 0, Improve]
    X[dat_spp[FeasChange < 0,SiteRef]] <- 0 - dat_spp[FeasChange < 0, Decline]
    X[X2<0.5] <- NA # remove cells where the feasibility expansion is less than 0.5 (X2 is from the mean feasibilty panel)
    
    breakpoints <- c(seq(-100, -50,10), seq(60, 100,10));length(breakpoints)
    labels <- c("Decline", "Improve")
    ColScheme <- c(brewer.pal(11,"RdBu")[c(1:4)], "grey90", "grey90", brewer.pal(11,"RdBu")[c(8:11)]); length(ColScheme)
    
    par(plt = c(0.6, 0.95, 0.25, 1), xpd = TRUE, new = TRUE)
    terra::image(X,xlab = NA,ylab = NA,bty = "n", xaxt="n", yaxt="n", col=ColScheme, breaks=breakpoints, asp = 1)
    terra::plot(outline, add=T, border="black",col = NA, lwd=0.4)
    
    xl <- 1600000; yb <- 1000000; xr <- 1700000; yt <- 1700000
    rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme)
    text(rep(xr+10000,length(labels)),seq(yb,yt,(yt-yb)/(15-1))[c(3,9)],labels,pos=4,cex=0.7,font=0.7, srt=90)
    text(rep(xr-20000,length(labels)),seq(yb,yt,(yt-yb)/(15-1))[c(1,8,15)],c("100%", "50%", "100%"),pos=4,cex=0.7,font=1)
    text(xl-30000, mean(c(yb,yt))-30000, paste("Ensemble agreement\n(% of GCMs)", sep=""), srt=90, pos=3, cex=0.75, font=2)
    par(xpd=F)
  }
  
  if(save_png){
    dev.off()
  } #saves to the main ccissr folder
} 

#1) Adding terra:: explicitly before every image() and plot(outline,) fixed the error!

### Make plots: -----------
bc_ol <- vect("data-raw/data_tables/bc_outline.gpkg")

plot_2panel(con, spp = "Pl", edatope = "C4", period = "2021_2040", bgc_template = bgc_template, outline = bc_ol, three_panel = FALSE, save_png = TRUE)







#2) Plot the reference period maps separately -------------------

plot_historical_suit_map <- function(dbCon, spp, edatope, bgc_template, outline, save_png = TRUE, out_dir = file.path("tree_summaries_code", "Figures"), width = 6, height = 5, res = 300) 
  {
  
  ## ===============================
  ## Output device
  if (save_png) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    out_file <- file.path(out_dir, paste("Historic_Suit_Map", spp, edatope, "png", sep = "."))
    png(filename = out_file, type = "cairo", units = "in", width = width, height = height, res = res)
    on.exit(dev.off(), add = TRUE)
  }
  
  ## ===============================
  ## Query & prep (NO FuturePeriod filter)
  dat_spp <- DBI::dbGetQuery(
    dbCon,
    sprintf(
      "SELECT SiteRef, Curr
       FROM cciss_res
       WHERE Spp = '%s'
         AND Edatope = '%s'",
      spp, edatope
    )
  ) |> data.table::as.data.table()
  
  # Classify Curr like the original workflow
  dat_spp[, Curr := as.integer(round(Curr))]
  dat_spp[is.na(Curr) | Curr > 3.5, Curr := 5]  #5 = "not mapped"/white
  
  ## ===============================
  ## Build raster
  X <- terra::rast(bgc_template$bgc_rast)
  terra::values(X) <- NA
  
  X[dat_spp$SiteRef] <- dat_spp$Curr
  
  ## ===============================
  ## Plot
  par(mar = c(0.5, 0.5, 0.5, 0.5), bty = "n", pty = "m")
  
  breakseq  <- c(0.5, 1.5, 2.5, 3.5, 5)
  ColScheme <- c("darkgreen", "dodgerblue1", "gold2", "white")
  
  terra::image(
    X,
    xlab = NA, ylab = NA,
    xaxt = "n", yaxt = "n",
    col = ColScheme,
    breaks = breakseq,
    asp = 1
  )
  
  terra::plot(outline, add = TRUE, border = "black", col = NA, lwd = 0.4)
  
  legend(
    "topright",
    legend = c("E1 (high)", "E2 (moderate)", "E3 (low)"),
    fill = ColScheme[1:3],
    bty = "n",
    cex = 1.0,
    title = "Historical suitability",
    inset = c(0.02, 0.02)
  )
  
  invisible(NULL)
}

## Plots: --------------
plot_historical_suit_map(dbCon = con, spp = "Pl", edatope = "D6", bgc_template = bgc_template, outline = bc_ol)













# 3) Modifications to the future-period map + boxplot ----------------------

plot_future_plus_boxplot <- function(
    dbCon,
    spp,
    edatope,
    period,
    bgc_template,
    outline,
    save_png = TRUE
) {
  
  ## ===============================
  ## Lookups and constants
  zoneScheme <- c(
    PP = "#ea7200", MH = "#6f2997", SBS = "#2f7bd2", ESSF = "#ae38b8",
    CWH = "#488612", BWBS = "#4f54cf", CWF = "#7577e7", IGF = "#77a2eb",
    CMX = "#71d29e", BG = "#dd1320", IDF = "#e5d521", MS = "#e44ebc",
    SWB = "#a1dbde", CRF = "#af3a13", WJP = "#73330e", ICH = "#1fec26",
    CDF = "#edf418", JPW = "#96b3a5", CMA = "#eae1ee", SBPS = "#6edde9",
    IMA = "#e3f1fa", GBD = "#4d433f", OW = "#582511", BAFA = "#eee4f1",
    MMM = "#FF00FF", MHRF = "#2612dc", MGP = "#f0aeab", FG = "#92696c",
    SGP = "#cca261", GO = "#f0a325", SBAP = "#51d5a7", IWF = "#d44273",
    BSJP = "#424160", MSSD = "#dac370", MDCH = "#2d0cd4", CVG = "#c9edd3",
    SAS = "#92b1b6", CCH = "#7e22ca"
  )
  
  spps.lookup   <- data.table::copy(ccissr::T1)
  edatope.names <- c("Poor-subxeric", "Medium-mesic", "Rich-hygric")
  edatopes      <- c("B2", "C4", "D6")
  
  ## ===============================
  ## Output device
  if (save_png) {
    
    out_dir <- file.path("tree_summaries_code", "Figures")
    
    out_file <- file.path(
      out_dir,
      paste("Future_Map_Box", spp, edatope, period, "png", sep = ".")
    )
    
    png(
      filename = out_file,
      type = "cairo",
      units = "in",
      width = 7,
      height = 4.5,
      res = 300
    )
    
    on.exit(dev.off(), add = TRUE)
  }
  
  ## ===============================
  ## Layout: text | map
  ##         box  | map
  layout(
    matrix(c(1, 3,
             2, 3),
           nrow = 2, byrow = TRUE),
    widths  = c(0.35, 0.65),
    heights = c(0.28, 0.72)
  )
  
  ## ===============================
  ## Data query and prep
  dat_spp <- dbGetQuery(
    dbCon,
    sprintf(
      "
      SELECT *
      FROM cciss_res
      WHERE Spp = '%s'
        AND Edatope = '%s'
        AND FuturePeriod = '%s'
      ",
      spp, edatope, period
    )
  ) |> data.table::as.data.table()
  
  dat_spp[, FeasChange := Curr - Newsuit] #Curr is the suitability from the reference period
  dat_spp[, Curr := as.integer(round(Curr))]
  dat_spp[is.na(Curr) | Curr > 3.5, Curr := 5]
  
  ## ===============================
  ## Build rasters for future change map
  X  <- terra::rast(bgc_template$bgc_rast)
  X2 <- terra::rast(X)
  X3 <- terra::rast(X)
  
  terra::values(X)  <- NA
  terra::values(X2) <- NA
  terra::values(X3) <- NA
  
  X[ dat_spp[Curr < 3.5 | Newsuit < 3.5, SiteRef] ] <-
    dat_spp[Curr < 3.5 | Newsuit < 3.5, FeasChange]
  
  X2[ dat_spp[Curr > 3.5 & Newsuit < 3.5, SiteRef] ] <-
    dat_spp[Curr > 3.5 & Newsuit < 3.5, FeasChange]
  
  X3[ dat_spp[Curr < 4 & Newsuit > 3.5, SiteRef] ] <- 1
  
  ## ===============================
  ## PANEL 1: Text (top-left)
  par(mar = c(0, 0.4, 0, 0))
  plot.new()
  
  Common <- spps.lookup[TreeCode == spp, EnglishName]
  
  mtext(
    bquote(bold(.(spp)) ~ "-" ~ .(Common)),
    side = 3, line = -1.5, adj = 0, cex = 1
  )
  
  if (edatope %in% edatopes) {
    mtext(
      paste(
        "Site type:", edatope,
        "(", edatope.names[edatope == edatopes], ")"
      ),
      side = 3, line = -3.2, adj = 0, cex = 0.85
    )
  }
  
  mtext(
    paste("Time period:", period),
    side = 3, line = -4.8, adj = 0, cex = 0.85
  )
  
  ## ===============================
  ## PANEL 2: Boxplot (bottom-left)
  par(mar = c(5, 2, 1, 0))
  
  bgc_mapped <- as.data.frame(bgc_template$bgc_rast, cells = TRUE)
  data.table::setDT(bgc_mapped)
  bgc_mapped[bgc_template$ids, bgc := i.bgc, on = "bgc_id"]
  bgc_mapped[, zone := regmatches(bgc, regexpr("^[A-Z]+", bgc))]
  
  dat_spp[bgc_mapped, zone := i.zone, on = c(SiteRef = "cell")]
  
  z <- boxplot(FeasChange ~ zone, data = dat_spp, plot = FALSE)
  
  for (i in seq_along(z$names)) {
    vals <- dat_spp[zone == z$names[i], FeasChange]
    z$stats[c(1, 5), i] <- quantile(vals, c(0.05, 0.95), na.rm = TRUE)
  }
  
  bxp(
    z,
    boxfill = zoneScheme[z$names],
    outline = FALSE,
    ylim = c(-3, 3),
    las = 2,
    cex.axis = 0.7
  )
  
  abline(h = 0, col = "grey50", lwd = 2)
  mtext("Mean change in suitability by BGC zone", side = 3, cex = 0.75, font = 2)
  
  ## ===============================
  ## PANEL 3: Future change map (right)
  par(mar = c(1, 0.2, 1, 1), bty = "n")
  
  breakpoints <- seq(-3, 3, 0.5)
  
  ColScheme <- c(
    "black",
    RColorBrewer::brewer.pal(11, "RdBu")[1:4],
    "grey90", "grey90",
    RColorBrewer::brewer.pal(11, "RdBu")[7:11]
  )
  
  ColScheme2 <- c(
    brewer.pal(11,"RdBu")[1:4],
    "grey90","grey90",
    colorRampPalette(c("white","khaki1","gold"))(6)
  )
  
  terra::image(X, col = ColScheme, breaks = breakpoints, asp = 1, xaxt = "n", yaxt = "n", xlab = NA, ylab = NA)
  terra::image(X2, add = TRUE, col = ColScheme2, breaks = breakpoints)
  terra::image(X3, add = TRUE, col = 1)
  
  terra::plot(outline, add = TRUE, col = NA, lwd = 0.4)
  
  ## ===============================
  ## LEGEND (top-right) — rectangles + text (base graphics)
  draw_future_legend <- function(ColScheme, ColScheme2,
                                 title = "Mean change in suitability") {
    
    usr <- par("usr")
    dx  <- usr[2] - usr[1]
    dy  <- usr[4] - usr[3]
    
    # Increase wf/hf to scale up; nudge x0f/y0f to keep it tucked in the corner.
    wf <- 0.42
    hf <- 0.65
    x0f <- 0.67
    y0f <- 0.40
    
    x0 <- usr[1] + x0f * dx
    y0 <- usr[3] + y0f * dy
    w  <- wf  * dx
    h  <- hf  * dy
    x1 <- x0 + w
    y1 <- y0 + h
    
    op <- par(xpd = NA)
    on.exit(par(op), add = TRUE)
    
    # Whole legend box 
    rect(x0, y0, x1, y1, col = NA, border = NA)
    
    # Internal padding
    px <- 0.05 * w
    py <- 0.05 * h
    
    # Column widths (fractions inside legend)
    w_title   <- 0.18 * w
    w_exp_bar <- 0.10 * w
    w_gap     <- 0.035 * w
    w_mean    <- 0.12 * w
    w_text    <- w - (2*px + w_title + w_exp_bar + w_gap + w_mean)
    
    # Vertical region for the bars (leave space below for Loss)
    bar_top    <- y1 - py
    bar_bottom <- y0 + 0.27 * h
    bar_h      <- bar_top - bar_bottom
    
    # Column x positions
    x_title_center <- x0 + px + 0.5*w_title
    x_exp0         <- x0 + px + w_title
    x_mean0        <- x_exp0 + w_exp_bar + w_gap
    x_text0        <- x_mean0 + w_mean
    
    # Slightly larger text overall
    cex_title <- 0.95
    cex_lab   <- 0.95
    cex_loss  <- 0.95
    
    # Title (vertical)
    text(x_title_center, (bar_top + bar_bottom)/2, title,
         srt = 90, font = 2, cex = cex_title)
    
    # --- Expansion column: gradient on top, white label box below
    # In the original: label box is BELOW the yellow gradient.
    exp_label_frac <- 0.36              # portion of bar height used by label box
    exp_lbl_h      <- exp_label_frac * bar_h
    exp_grad_h     <- bar_h - exp_lbl_h
    
    # Gradient (top part)
    exp_cols <- rev(tail(ColScheme2, 6))  # reverse for darkest at top
    bin_h_e  <- exp_grad_h / 6
    for (i in 1:6) {
      yb0 <- bar_top - i*bin_h_e
      yb1 <- yb0 + bin_h_e
      rect(x_exp0, yb0, x_exp0 + w_exp_bar, yb1,
           col = exp_cols[i], border = "black", lwd = 1)
    }
    
    # White label box (bottom part)
    rect(x_exp0, bar_bottom, x_exp0 + w_exp_bar, bar_bottom + exp_lbl_h,
         col = "white", border = "black", lwd = 1)
    text(x_exp0 + 0.5*w_exp_bar, bar_bottom + 0.5*exp_lbl_h, "Expansion",
         srt = 90, cex = cex_lab)
    
    # --- Mean-change bar 
    mean_cols <- c(
      ColScheme[length(ColScheme)],     # +3
      ColScheme[length(ColScheme)-1],   # +2
      ColScheme[length(ColScheme)-2],   # +1
      "grey90",                         # no change
      ColScheme[4],                     # -1
      ColScheme[3],                     # -2
      "black"                           # -3
    )
    bin_h <- bar_h / 7
    for (i in 1:7) {
      yb0 <- bar_top - i*bin_h
      yb1 <- yb0 + bin_h
      rect(x_mean0, yb0, x_mean0 + w_mean, yb1,
           col = mean_cols[i], border = "black", lwd = 1)
    }
    
    # Labels for mean-change bins
    labs <- c("+3", "+2", "+1", "no change", "-1", "-2", "-3")
    for (i in 1:7) {
      ymid <- bar_top - (i - 0.5)*bin_h
      text(x_text0 + 0.01*w, ymid, labs[i], adj = 0, cex = cex_lab)
    }
    
    # --- Loss swatch under bars
    gap_below_bars <- 0.02 * h  # vertical gap between bars and Loss
    loss_h <- bin_h   # match height of the 7-bin rectangles
    loss_yc <- (bar_bottom - gap_below_bars) - loss_h/2  # center just below bars
    
    rect(x_mean0, loss_yc - loss_h/2, x_mean0 + w_mean, loss_yc + loss_h/2,
         col = "black", border = "black", lwd = 1)
    text(x_text0 + 0.01*w, loss_yc, "Loss", adj = 0, cex = cex_loss)
  }
  
  draw_future_legend(ColScheme = ColScheme, ColScheme2 = ColScheme2)
  
  invisible(NULL)
}





## Plots: ------------
plot_future_plus_boxplot(dbCon = con,
  spp = "Pl",
  edatope = "B2",
  period = "2081_2100", # quick grab: 2021_2040, 2041_2060, 2061_2080, 2081_2100
  bgc_template = bgc_template, outline = bc_ol,
  save_png = TRUE 
)
