
##======================================================================================
## CCISS Publication Scripts
## Step 4g - Figures - Manuscript plot of Suitability and species persistence
##======================================================================================

# Colin Mahony
# c_mahony@alumni.ubc.ca
# 778-288-4008
# July 21, 2019




################################
## Manuscript figure
###################################

edatope=edatopes[2]
rcp=rcps[1]
proj.year=proj.years[2]

for(proj.year in proj.years){
  for(edatope in edatopes){
    SuitRichness.ref <- get(paste("SuitRichness.ref", edatope, sep="."))
    SuitRichness.proj <- get(paste("SuitRichness", rcp, proj.year, edatope, sep="."))
    SuitRichnessChange <- get(paste("SuitRichnessChange", rcp, proj.year, edatope, sep="."))
    SuitRichnessChangePct <- SuitRichness.proj/SuitRichness.ref
    SuitRichnessChangePct[!is.finite(SuitRichnessChangePct)] <- NA
    SppRichness.ref <- get(paste("SppRichness.ref", edatope, sep="."))
    SppRichness.proj <- get(paste("SppRichness", rcp, proj.year, edatope, sep="."))
    SppRichnessChange <- get(paste("SppRichnessChange", rcp, proj.year, edatope, sep="."))
    SppRichnessChangePct <- SppRichness.proj/SppRichness.ref
    SppRichnessChangePct[!is.finite(SppRichnessChangePct)] <- NA
    SuitTurnover <- get(paste("SuitTurnover", rcp, proj.year, edatope, sep="."))
    SppTurnover <- get(paste("SppTurnover", rcp, proj.year, edatope, sep="."))
    SuitTurnover[which(SuitTurnover==99)] <- NA
    SppTurnover[which(SppTurnover==99)] <- NA
    SuitPersistence <- get(paste("SuitPersistence", rcp, proj.year, edatope, sep="."))
    SppPersistence <- get(paste("SppPersistence", rcp, proj.year, edatope, sep="."))
    
    metric <- "SuitPersistence"
    
    # x11(width=6.5, height=5, pointsize=8)
    png(filename=paste("results\\Manu_Persistence\\CCISS_manu_", metric, edatope, rcp, proj.year,"png",sep="."), type="cairo", units="in", width=6.5, height=4.55, pointsize=8, res=400)
    # pdf(file=paste("Results\\CCISS.Fig7", metric, edatope, rcp, proj.year,"pdf",sep="."), width=6.5, height=4.55, pointsize=8)
    
    ylim <- c(0,1.4)
    y <- get(metric)
    # y[y<2^(ylim[1])] <- 2^(ylim[1])
    # y <- log2(y)
    
    par(mar=c(0,0,0,0))
    plot(0, xlim=c(0,1), ylim=c(0,1), col="white", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
    # box()
    
    #map
    par(mar=c(0,0,0,0), plt = c(0.05, 0.85, 0.025, 0.999), new = TRUE)
    # SuitTurnover <- get(paste("SuitTurnover", rcp, proj.year, edatope, sep="."))
    # SppTurnover <- get(paste("SppTurnover", rcp, proj.year, edatope, sep="."))
    values(X) <- y[plotOrder]
    
    breakpoints <- seq(0,2,0.2); length(breakpoints)
    labels <- c("0%", "50%", "No change", "150%", "200%")
    ColScheme <- c(brewer.pal(11,"RdBu")[1:4], rep("gray90", 2), brewer.pal(11,"RdBu")[8:11]); length(ColScheme)
    # ColScheme <- brewer.pal(11,"RdBu")[-6]; length(ColScheme)
    
    plot(bdy.bc, border="black", lwd=0.4)
    image(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks=breakpoints, maxpixels= ncell(X))
    mtext("(a)", side=3, line=-2.5, adj=0.0, cex=1, font=2)
    par(xpd=T)
    xl <- 250000; yb <- 1000000; xr <- 350000; yt <- 1500000
    # xl <- 2025000; yb <- 500000; xr <- 2100000; yt <- 950000
    rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme)
    text(rep(xr,length(labels)),seq(yb,yt,(yt-yb)/(length(labels)-1)),labels,pos=4,cex=1,font=1)
    # text(rep(xl,2),c(yb,yt)+c(-10000, 10000),c("No change", "Full turnover"),pos=c(1,3),cex=1,font=1)
    text(xl-40000, mean(c(yb,yt))-30000, paste("Feasibility persistence\n(", proj.year.name[which(proj.years==proj.year)], ", ", c("RCP4.5", "RCP8.5")[which(rcps==rcp)], ")", sep=""), srt=90, pos=3, cex=1, font=2)
    # rect(xl,  yt+20000,  xr,  yt+60000,  col=ColScheme[length(ColScheme)])
    # text(xr,  yt+40000,  bquote(">"*.(breakseq[3])*sigma),pos=4,cex=1,font=1)  
    par(xpd=F)
    mtext(paste(edatope.name[which(edatopes==edatope)], " sites (", edatope, ")", sep=""), side=3, line=-1.5, adj=0.35, cex=1, font=2)
    
    ##==============================
    ## Callout box for Location 1
    bgc.select <- "MSdk"
    q.select <- 0.91 # 0.90 is good
    focal.bgc <- which(BGC==bgc.select)
    pt <- focal.bgc[which(SuitTurnover[focal.bgc]==min(SuitTurnover[focal.bgc][SuitTurnover[focal.bgc]>quantile(SuitTurnover[focal.bgc], q.select, na.rm=T)], na.rm=T))[1]] #select the grid point
    pt <- pt[1]
    # pt <- 96923
    lines(c(pts.df[pt,1], pts.df[pt,1]+200000), c(pts.df[pt,2], pts.df[pt,2]+0), lwd=1.5)
    points(pts.df[pt,], pch=21, cex=1.75, lwd=1.5, bg=ColScheme[cut(y[pt],breaks=breakpoints)])
    
    # assemble the species mix for the reference period
    comm.ref.pt <- rep(NA, length(spps))
    for(spp in spps){
      Suit <- read.csv(paste("outputs\\Suit.ref", grid, spp, edatope, "csv", sep="."))[pt,1]
      Suit[is.na(Suit)] <- 5  #XXX note this is different from the equivalent line for the other time periods.
      Suit <- 1-(Suit-1)/4
      comm.ref.pt[which(spps==spp)] <- Suit
    }
    names(comm.ref.pt) <- spps
    #convert fractional suitability to standard suitability
    comm.ref.pt <- round((1-comm.ref.pt)*4+1)
    comm.ref.pt[comm.ref.pt==4] <- 5
    
    # assemble the projected community
    for(GCM in GCMs){
      comm.proj <- rep(NA, length(spps))
      for(spp in spps){
        Suit <- read.csv(paste("outputs\\Suit", grid, GCM, rcp, proj.year, spp, edatope, "csv", sep="."))[pt,1]
        Suit[Suit==5] <- 5
        Suit <- 1-(Suit-1)/4
        comm.proj[which(spps==spp)] <- Suit
      }
      names(comm.proj) <- spps
      if(GCM==GCMs[1]) comm.proj.pt <- comm.proj else comm.proj.pt <- rbind(comm.proj.pt,comm.proj)
      # print(GCM)
    }
    comm.proj.pt.ensMean <- apply(comm.proj.pt, 2, FUN=mean, na.rm=T)
    comm.proj.pt.ensMean <- round((1-comm.proj.pt.ensMean)*4+1)
    comm.proj.pt.ensMean[comm.proj.pt.ensMean==4] <- 5
    
    std.ref <- paste(paste(names(comm.ref.pt)[which(comm.ref.pt==1)],collapse=""), "(", 
                     paste(names(comm.ref.pt)[which(comm.ref.pt==2)],collapse=""), ")((", 
                     paste(names(comm.ref.pt)[which(comm.ref.pt==3)],collapse=""), "))", sep="")
    std.proj <- paste(paste(names(comm.proj.pt.ensMean)[which(comm.proj.pt.ensMean==1)],collapse=""), "(", 
                      paste(names(comm.proj.pt.ensMean)[which(comm.proj.pt.ensMean==2)],collapse=""), ")((", 
                      paste(names(comm.proj.pt.ensMean)[which(comm.proj.pt.ensMean==3)],collapse=""), "))", sep="")
    
    BGC.list <- vector()
    for(GCM in GCMs){
      BGC.list[which(GCMs==GCM)] <- get(paste("BGC.pred", GCM, rcp, proj.year, sep="."))[pt]
      # print(BGC.pred[pt])
      # print(GCM)
    }
    BGC.list 
    
    par(mar=c(0,0,0,0), plt = c(0.825, 0.995, 0.01, 0.445), new = TRUE, mgp=c(2,0.1,0))
    plot(0, xlim=c(0,1), ylim=c(0,1), col="white", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
    rect(-9,-9,9,9, col="white", lwd=0)
    mtext("(d)", side=3, line=-1, adj=-0.2, cex=1, font=2)
    box()
    
    # par(mar=c(0,0,0,0), plt = c(0.6, 0.99, 0.01, 0.175), new = TRUE, mgp=c(2,0.1,0))
    # plot(0, xlim=c(0,1), ylim=c(0,1), col="white", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
    pos1 <- 1
    text(-0.025, pos1-0, paste("Historical bioclimate"), pos=4, cex=0.8, font=2, offset=0.1)
    text(0, pos1-0.05, bgc.select, pos=4, cex=0.8, offset=0.1)
    text(-0.025, pos1-.125, paste("Projected bioclimates"), pos=4, cex=0.8, font=2, offset=0.1)
    pos2 <- 0.44
    text(-.025, pos2-0, "Historical feasibilities", pos=4, cex=0.8, font=2, offset=0.1)
    text(0.0, pos2-0.05, std.ref, pos=4, cex=0.8, offset=0.1)
    text(-.025, pos2-0.125, paste("Projected feasibilities", sep=""), pos=4, cex=0.8, font=2, offset=0.1)
    for(i in 1:length(unique(BGC.list))){
      BGC.i <- names(rev(sort(table(BGC.list))))[i]
      comm.proj.pt.i <- comm.proj.pt[which(BGC.list==BGC.i)[1],]
      std.proj.i <- paste(paste(names(comm.proj.pt.i)[which(comm.proj.pt.i==1)],collapse=""), "(", 
                          paste(names(comm.proj.pt.i)[which(comm.proj.pt.i==0.75)],collapse=""), ")((", 
                          paste(names(comm.proj.pt.i)[which(comm.proj.pt.i==0.5)],collapse=""), "))", sep="")
      text(0.0, pos2-0.125-0.05*i, paste(BGC.i, ": ", std.proj.i, sep=""), pos=4, cex=0.75, offset=0.1)
    }
    
    text(-.025, 0, paste("Mean persistence: ", round(SuitPersistence[pt]*100), "%", sep=""), pos=4, cex=0.8, offset=0.1, font=2)
    
    par(mar=c(3,0,1,0), plt = c(0.855, 0.99, 0.305, 0.36), new = TRUE, mgp=c(0.5,0.2,0))
    bar.bgc <- names(rev(sort(table(BGC.list))))
    bar.zone <- rep(NA, length(bar.bgc))
    for(i in BGCcolors$classification){ bar.zone[grep(i,bar.bgc)] <- i }
    barplot(rev(sort(table(BGC.list))), col=as.character(BGCcolors$colour[match(bar.zone, BGCcolors$classification)]), horiz=F, las=2, ylab=list("# GCMs", cex=0.7), yaxt="n", cex.names=0.7)
    par(mgp=c(1,0,0))
    axis(2, at=seq(0,12,4), labels=seq(0,12,4), las=1, tck=0, cex.axis=0.8, lwd=0)
    
    # spp.ref <- c(names(comm.ref.pt)[which(comm.ref.pt==1)], names(comm.ref.pt)[which(comm.ref.pt==2)], names(comm.ref.pt)[which(comm.ref.pt==3)])
    # boxplot(0, col="white", xlim=c(1,length(spp.ref)), ylim=c(0,1))
    # for(spp in spp.ref){
    #   y <- comm.proj.pt[,which(names(comm.proj.pt)==spp)]
    #   boxplot(y, add=T, at=which(spp.ref==spp), range=0, names=spp)
    # }
    # axis(1, at=1:length(spp.ref), labels = spp.ref, tck=0)
    
    
    ##===========================
    ## Summary by zone
    par(mar=c(4.5,2,0.1,0.1), plt = c(0.1, 0.325, 0.125, 0.325), new = TRUE, mgp=c(2.5,0.25,0))
    xlim=c(1, length(levels(droplevels(zone[-grep("BAFA|CMA|IMA",zone)]))))
    z <- boxplot(y~zone, ylim=ylim, ylab="", vertical = TRUE, plot=F)
    for(i in 1:length(levels(zone))){ 
      temp <- y[which(zone==levels(zone)[i])]
      z$stats[c(1,5), i] <- quantile(temp[!is.na(temp)],c(0.05, 0.95))
    }
    bxp(z, ylim=ylim, xlim=xlim, xaxt="n", yaxt="n", xaxs="i", ylab="", pch=0,outline=FALSE)
    lines(c(-99,99), c(1,1), lwd=2, col="darkgrey")
    bxp(z, add=T, boxfill = as.character(BGCcolors$colour[match(levels(droplevels(zone[-grep("BAFA|CMA|IMA",zone)])), BGCcolors$classification)]), xaxt="n", yaxt="n",xaxs="i", 
        ylab=paste("Feasibility persistence\n(", proj.year.name[which(proj.years==proj.year)], ", ", c("RCP4.5", "RCP8.5")[which(rcps==rcp)], ")", sep=""),  pch=0,outline=FALSE)
    axis(1, at=1:length(levels(zone)), levels(zone), tick=F, las=2, cex=0.8)
    axis(2, seq(ylim[1],ylim[2],.25), paste(seq(ylim[1],ylim[2],.25)*100, "%", sep=""), tick=F, las=2)
    # lines(range(grep("CWH|MH|CDF", levels(zone)))+c(-0.6,0.6), if(scenario=="6variable") rep(2.7,2) else rep(3.2,2), col=alpha("red", 0.5), lwd=3, lty=1)
    # lines(range(grep("ESSF|MS|IDF|PP|BG|ICH|SBPS|SBS|BWBS|SWB", levels(zone)))+c(-0.6,0.6), if(scenario=="6variable") rep(1.5,2) else rep(1.6,2), col=alpha("red", 0.5), lwd=3, lty=1)
    axis(1, at=c(mean(grep("CWH|MH|CDF", levels(zone))), mean(grep("ESSF|MS|IDF|PP|BG|ICH|SBPS|SBS|BWBS|SWB", levels(zone))),mean(grep("CMA|IMA|BAFA", levels(zone)))), 
         labels = c("Coast", "Interior", "Alpine"), tick=F, line=2.95 )
    par(xpd=NA)
    x1 <- range(grep("CWH|MH|CDF", levels(zone)))+c(-0.4,0.4)
    x2 <- range(grep("ESSF|MS|IDF|PP|BG|ICH|SBPS|SBS|BWBS|SWB", levels(zone)))+c(-0.4,0.4)
    # x3 <- range(grep("CMA|IMA|BAFA", levels(zone)))+c(-0.4,0.4)
    bracketpos=-0.7
    lines(rep(x1, each=2), bracketpos*c(0.95, 1,1,0.95))
    lines(rep(x2, each=2), bracketpos*c(0.95, 1,1,0.95))
    # lines(rep(x3, each=2), bracketpos*c(0.95, 1,1,0.95))
    par(xpd=F)  
    mtext("(b)", side=3, line=0.25, adj=-0.15, cex=1, font=2)
    
    
    #===============================================================================
    # Suitability Persistence relative to temperature change (predicted baseline)
    #===============================================================================
    library(msir)
    ColScheme=c(2,1,3)
    
    plt=c(0.665, 0.995, 0.525, 0.99)
    ylim=c(0,1.1)
    par(mfrow=c(1,1), mar=c(3.25,3.25,0.1,0.1), mgp=c(1.5,0.25,0), plt = plt, new = TRUE)
    plot(0, xlim=c(0,7.9), ylim=ylim, yaxs="i", xaxs="i", col="white", xaxt="n", yaxt="n", 
         xlab=bquote(BC~mean~temperature~change~relative~to~"1961-90"~"("*degree*C*")"), 
         ylab="")
    axis(1, at=0:8, labels = 0:8, tck=0)
    axis(2, at=seq(0,1.2,0.2), labels = paste(seq(0,1.2,0.2)*100, "%", sep=""), las=2, tck=0)
    par(mgp=c(2,0.25,0), plt = plt, new = TRUE)
    title(ylab="BC mean feasibility persistence             ")
    mtext("(c)", side=3, line=-0.75, adj= -0.1, cex=1, font=2)
    lines(c(-99,99), c(1,1), lwd=2, col="darkgrey")
    
    for(edatope.temp in edatopes){
      x <- c(0,MAT.change)
      y <- c(1,get(paste(metric, "mean.all", edatope.temp, sep=".")))
      l <- loess.sd(y~x, span=1, nsigma = 1)
      if(edatope.temp!=edatope) polygon(c(l$x, rev(l$x)), c(l$upper, rev(l$lower)), col=alpha(ColScheme[which(edatopes==edatope.temp)], 0.5), border=NA)
      l <- loess(y~x)
      par(xpd=T)
      text(max(x)-0.1, predict(l, max(x)), edatope.temp, pos=4, font= if(edatope.temp==edatope) 2 else 1, col=ColScheme[which(edatopes==edatope.temp)], cex= if(edatope.temp==edatope) 1.1 else 1)
      par(xpd=F)
    }
    
    temp.suit <- rep(NA, length(GCMs))
    temp.spp <- rep(NA, length(GCMs))
    for(rcp.temp in rcps){
      for(proj.year.temp in proj.years){
        temp.suit <- rep(NA, length(GCMs))
        temp.spp <- rep(NA, length(GCMs))
        
        x <- get(paste("MAT.change", rcp.temp, proj.year.temp, sep="."))
        y <- get(paste(metric, "mean", rcp.temp, proj.year.temp, edatope, sep=".")) 
        
        points(x,y, pch=c(21,22)[which(rcps==rcp.temp)], bg=c("black", "dodgerblue", "yellow")[which(proj.years==proj.year.temp)])
        # points(x,y2, col="red", pch=16)
        # points(0,mean(BGC.pred.ref!=BGC), pch=16, cex=1.3)
        # points(0,mean(zone.pred.ref!=zone), col="red", pch=16, cex=1.3)
      }
    }
    
    for(hist.year in hist.years[c(4)]){
      x <- get(paste("MAT.change", hist.year, sep="."))
      # y <- get(paste("SuitPersistence.mean", hist.year, edatope, sep="."))
      y <- get(paste(metric, "mean", hist.year, edatope, sep="."))
      points(x,y, cex=1.4, pch=2)
      text(x,y-0.01, paste(hist.year.name[which(hist.years==hist.year)], ""), pos=1, cex=0.8, font=2)
      # print(y)
    }
    
    legend("bottomleft", legend=c("RCP4.5", "RCP8.5", "2011-2040", "2041-2070", "2071-2100"), pch=c(21,22, NA,NA,NA), 
           pt.bg=c("gray", "gray", NA,NA,NA), pt.cex=1.5, fill=c(NA, NA, "black", "dodgerblue", "yellow"), border = c(F, F, T, T, T), bty="n")
    
    # boxplot for focal period
    x <- c(0,MAT.change)
    x.focal <- MAT.change[which(seq.rcp==rcp & seq.proj.year==proj.year)]
    boxplot(x.focal, add=T, horizontal=TRUE, axes=FALSE, range=0, at=ylim[2]-0.05, boxwex = 0.12)
    text(max(x.focal), ylim[2]-0.05, paste(rcp.name[which(rcps==rcp)], ", ", proj.year.name[which(proj.years==proj.year)], sep=""), pos=4, cex=0.8)
    
    dev.off()
    print(edatope)
  }
  print(proj.year)
}





################################
## Supplemental figure of edatopes and time periods
###################################

rcp=rcps[1]

metric <- "SppPersistence"

# x11(width=6.5, height=5, pointsize=8)
png(filename=paste("results\\Manu_Persistence\\CCISS_manu_", metric,"png",sep="."), type="cairo", units="in", width=6.5, height=8.5, pointsize=11, res=400)
# pdf(file=paste("results\\CCISS_SummaryByBGC_", metric,".pdf",sep=""),  width=7.5, height=5.625, pointsize=15)
par(mar=c(0.1,0.1, 0.1,0.1), mgp=c(2,0.25,0), mfrow=c(3,2))

for(edatope in edatopes){
  for(proj.year in proj.years[1:2]){
    SuitPersistence <- get(paste("SuitPersistence", rcp, proj.year, edatope, sep="."))
    SppPersistence <- get(paste("SppPersistence", rcp, proj.year, edatope, sep="."))

    y <- get(metric)
    values(X) <- y[plotOrder]
    
    breakpoints <- seq(0,2,0.2); length(breakpoints)
    labels <- c("0%", "50%", "No change", "150%", "200%")
    ColScheme <- c(brewer.pal(11,"RdBu")[1:4], rep("gray90", 2), brewer.pal(11,"RdBu")[8:11]); length(ColScheme)
    # ColScheme <- brewer.pal(11,"RdBu")[-6]; length(ColScheme)
    
    plot(bdy.bc, border="black", lwd=0.4)
    image(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks=breakpoints, maxpixels= ncell(X))
    par(xpd=T)
    xl <- 1550000; yb <- 1000000; xr <- 1650000; yt <- 1700000
    # xl <- 2025000; yb <- 500000; xr <- 2100000; yt <- 950000
    rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme)
    text(rep(xr,length(labels)),seq(yb,yt,(yt-yb)/(length(labels)-1)),labels,pos=4,cex=1,font=1)
    # text(rep(xl,2),c(yb,yt)+c(-10000, 10000),c("No change", "Full turnover"),pos=c(1,3),cex=1,font=1)
    text(xl-40000, mean(c(yb,yt))-30000, paste("Feasibility persistence\n(", proj.year.name[which(proj.years==proj.year)], ", ", c("RCP4.5", "RCP8.5")[which(rcps==rcp)], ")", sep=""), srt=90, pos=3, cex=1, font=2)
    # rect(xl,  yt+20000,  xr,  yt+60000,  col=ColScheme[length(ColScheme)])
    # text(xr,  yt+40000,  bquote(">"*.(breakseq[3])*sigma),pos=4,cex=1,font=1)  
    par(xpd=F)
    mtext(paste("(", if(proj.year==proj.years[1]) letters[c(1,3,5)[which(edatopes==edatope)]] else letters[c(2,4,6)[which(edatopes==edatope)]], ")\n", proj.year.name[which(proj.years==proj.year)], "\n", edatope.name[which(edatopes==edatope)], " (", edatope, ")", sep=""), side=1, line=-1.5, adj=0.01, cex=1, font=2)
    
    print(proj.year)
  }
  print(edatope)
}
dev.off()



################################
## basic map
###################################

rcp=rcps[1]

metric <- "SuitPersistence"

for(edatope in edatopes){
  for(proj.year in proj.years[1:2]){

        y <- get(paste(metric, rcp, proj.year, edatope, sep="."))
    values(X) <- y[plotOrder]
    
    breakpoints <- seq(0,2,0.2); length(breakpoints)
    labels <- c("0%", "50%", "No change", "150%", "200%")
    ColScheme <- c(brewer.pal(11,"RdBu")[1:4], rep("gray90", 2), brewer.pal(11,"RdBu")[8:11]); length(ColScheme)
    # ColScheme <- brewer.pal(11,"RdBu")[-6]; length(ColScheme)
    
    # x11(width=6.5, height=5, pointsize=8)
    png(filename=paste("results\\Manu_Persistence\\CCISS_", metric,"_", rcp,"_", proj.year,"_", edatope,".png",sep=""), type="cairo", units="in", width=6.5, height=5.5, pointsize=11, res=400)
    # pdf(file=paste("results\\CCISS_SummaryByBGC_", metric,".pdf",sep=""),  width=7.5, height=5.625, pointsize=15)
    par(mar=c(0.1,0.1, 0.1,0.1), mgp=c(2,0.25,0))
    
    plot(bdy.bc, border="black", lwd=0.4)
    image(X, add=T, xaxt="n", yaxt="n", col=ColScheme, breaks=breakpoints, maxpixels= ncell(X))
    par(xpd=T)
    xl <- 1550000; yb <- 1000000; xr <- 1650000; yt <- 1700000
    # xl <- 2025000; yb <- 500000; xr <- 2100000; yt <- 950000
    rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  col=ColScheme)
    text(rep(xr,length(labels)),seq(yb,yt,(yt-yb)/(length(labels)-1)),labels,pos=4,cex=1,font=1)
    # text(rep(xl,2),c(yb,yt)+c(-10000, 10000),c("No change", "Full turnover"),pos=c(1,3),cex=1,font=1)
    text(xl-40000, mean(c(yb,yt))-30000, paste("Combined species feasibility\nrelative to 1961-1990"), srt=90, pos=3, cex=1, font=2)
    # rect(xl,  yt+20000,  xr,  yt+60000,  col=ColScheme[length(ColScheme)])
    # text(xr,  yt+40000,  bquote(">"*.(breakseq[3])*sigma),pos=4,cex=1,font=1)  
    par(xpd=F)
    mtext(paste(proj.year.name[which(proj.years==proj.year)], " (", rcp.name[which(rcps==rcp)], ")\n", edatope.name[which(edatopes==edatope)], " (", edatope, ") sites", sep=""), side=1, line=-2.5, adj=0.05, cex=1, font=2)
    plot(bdy.bc, border="black", lwd=0.4, add=T)
    
    dev.off()
    
    print(proj.year)
  }
  print(edatope)
}


