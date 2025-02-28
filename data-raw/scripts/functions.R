library(scales)
library(EnvStats)
library(plotly)

analog_novelty <- function(clim.targets, clim.analogs, label.targets, label.analogs, vars,
                           clim.icvs = NULL, label.icvs = NULL, weight.icv = 0.5, sigma = TRUE,
                           analog.focal = NULL, threshold = 0.95, pcs = NULL, 
                           plotScree = FALSE, 
                           plot2d = FALSE, plot2d.pcs = cbind(c(1,2,3,4), c(2,3,4,5)), 
                           plot3d = FALSE, plot3d.pcs=c(1,2,3), biplot = TRUE){
  
  analogs <- if(is.null(analog.focal)) unique(label.targets) else analog.focal # list of analogs to loop through
  novelty <- rep(NA, length(label.targets)) # initiate a vector to store the sigma dissimilarities
  
  for(analog in analogs){ # loop through all of the analogs used to describe the target climates. 
    clim.analog <- clim.analogs[label.analogs==analog, ..vars]
    clim.target <- clim.targets[label.targets==analog, ..vars]
    if(!is.null(clim.icvs)) clim.icv <- clim.icvs[label.icvs==analog, ..vars]
    
    ## data cleaning
    clim.analog <- clim.analog[complete.cases(clim.analog)] # remove rows without data
    clim.analog <- clim.analog[, .SD, .SDcols = which(sapply(clim.analog, function(x) var(x, na.rm = TRUE) > 0))]  # Remove zero-variance columns
    clim.target <- clim.target[, .SD, .SDcols = names(clim.analog)]
    if(!is.null(clim.icvs)) clim.icv <- clim.icv[complete.cases(clim.icv)]
    if(!is.null(clim.icvs)) clim.icv <- clim.icv[, .SD, .SDcols = names(clim.analog)]
    
    ## scale the data to the variance of the analog, since this is what we will ultimately be measuring the M distance in. 
    clim.mean <- clim.analog[, lapply(.SD, mean, na.rm = TRUE)]
    clim.sd <- clim.analog[, lapply(.SD, sd, na.rm = TRUE)]
    clim.analog[, (names(clim.analog)) := lapply(names(clim.analog), function(col) {
      (get(col) - unlist(clim.mean)[col]) / unlist(clim.sd)[col]
    })]
    clim.target[, (names(clim.target)) := lapply(names(clim.target), function(col) {
      (get(col) - unlist(clim.mean)[col]) / unlist(clim.sd)[col]
    })]
    if(!is.null(clim.icvs)) clim.icv[, (names(clim.icv)) := lapply(names(clim.icv), function(col) {
      (get(col) - unlist(clim.icv[, lapply(.SD, mean, na.rm = TRUE)])[col]) / unlist(clim.sd)[col] # subtract mean of ICV to centre the ICV on zero. 
    })]
    
    ## PCA on pooled target and analog
    s <- sample(1:dim(clim.target)[1], dim(clim.analog)[1], replace = TRUE) # select a random sample of the target population to match the analog points. bootstrap if target population is smaller than analog points
    clim.target.sample <- clim.target[s,]
    pca <- prcomp(rbind(clim.analog, clim.target.sample), scale=FALSE)
    pcs.analog <- data.table(predict(pca, clim.analog))
    pcs.target <- data.table(predict(pca, clim.target))
    if(!is.null(clim.icvs)) pcs.icv <- data.table(predict(pca, clim.icv))
    
    if(is.null(pcs)){
      ## select number of pcs
      cumvar <- cumsum(pca$sdev^2 / sum(pca$sdev^2)) # vector of cumulative variance explained
      pcs <- which(cumvar >= threshold)[1]
      if(pcs<3) pcs <- 3
    }
    
    ## z-standardize the pcs to the variance of the analog. this is necessary for a metric that can be translated into sigma values. 
    weight.analog <- 1 - weight.icv
    pcs.mean.analog <- pcs.analog[, lapply(.SD, mean, na.rm = TRUE)]
    pcs.sd.analog <- pcs.analog[, lapply(.SD, sd, na.rm = TRUE)]
    if(!is.null(clim.icvs)) pcs.sd.icv <- pcs.icv[, lapply(.SD, sd, na.rm = TRUE)]
    if(!is.null(clim.icvs)) pcs.sd.combined <- weight.analog * pcs.sd.analog + weight.icv * pcs.sd.icv
    pcs.sd.use <- if(!is.null(clim.icvs)) pcs.sd.combined else pcs.sd.analog
    pcs.analog[, (names(pcs.analog)) := lapply(names(pcs.analog), function(col) {
      (get(col) - unlist(pcs.mean.analog)[col]) / unlist(pcs.sd.use)[col]
    })]
    pcs.target[, (names(pcs.target)) := lapply(names(pcs.target), function(col) {
      (get(col) - unlist(pcs.mean.analog)[col]) / unlist(pcs.sd.use)[col]
    })]
    if(!is.null(clim.icvs)) pcs.icv[, (names(pcs.icv)) := lapply(names(pcs.icv), function(col) {
      (get(col) - unlist(pcs.icv[, lapply(.SD, mean, na.rm = TRUE)])[col]) / unlist(pcs.sd.use)[col] # separately centering on the ICV mean becuase sometime the ICV is not centred on the centroid, and we want it to be. 
    })]
    
    ## create a combined covariance matrix for spatial variation and ICV
    cov.analog <- var(pcs.analog[, 1:pcs])
    cov.icv <- if (!is.null(clim.icvs)) var(pcs.icv[, 1:pcs]) else NULL
    if (!is.null(cov.icv)) {
      cov.combined <- weight.analog * cov.analog + weight.icv * cov.icv
    } else {
      cov.combined <- cov.analog
    }
    
    ## Mahalanobis distance and sigma dissimilarity
    md <- (mahalanobis(pcs.target[,1:pcs], rep(0, pcs), cov.combined))^0.5
    p <- pchi(md,pcs) # percentiles of the M distances on the chi distribution with degrees of freedom equaling the dimensionality of the distance measurement (PCs)
    q <- qchi(p,1) # values of the chi percentiles on a standard half-normal distribution (chi distribution with one degree of freedom)
    q[!is.finite(q)] <- 8 # set infinite values to 8 sigma (outside the decimal precision of pchi) 
    q[is.na(p)] <- NA # reset NA values as NA
    
    ## populate the novelty vector
    novelty[label.targets==analog] <- if(sigma) q else md
    
  } # end of the for-loop
  
  ## Plots for the final iteration of the for loop
  
  # Color Scheme for sigma novelty
  breakseq <- c(0,4,8)
  breakpoints <- c(seq(breakseq[1], breakseq[3], 0.01),199); length(breakpoints)
  ColScheme <- c(colorRampPalette(c("gray90", "gray50", "#FFF200", "#CD0000", "black"))(length(breakpoints)))
  
  ## Scree plot
  if(plotScree){
    par(mfrow=c(1,1), mar=c(3,3,1,1), mgp=c(1.75,0.25,0))
    a <- apply(predict(pca, clim.analog), 2, sd)
    b <- apply(predict(pca, clim.target), 2, sd)
    if(!is.null(clim.icvs)) c <- apply(predict(pca, clim.icv), 2, sd)
    diff <- abs(apply(predict(pca, clim.target), 2, mean) - apply(predict(pca, clim.analog), 2, mean))
    plot(0, xlim=c(1,length(a)), ylim=c(0,max(c(a,b, diff))*1.02), yaxs="i", col="white", tck=-0.005,
         xlab="Principal Component (PC)", ylab = "Standard Deviation")
    rect(pcs+0.5, -99, 99, 99, col = "grey95", lty=2)
    points(a, pch=21, bg="dodgerblue", cex=1.6)
    points(b, bg="grey", pch=21, cex=1.3)
    if(!is.null(clim.icvs)) points(c, bg="black", pch=21, cex=1)
    points(diff, col="black", pch=17, cex=1.3)
    text(pcs+0.5, max(c(a,b, diff)), paste0("Truncation at ", pcs, " PCs"), pos=4)
    s <- if(!is.null(clim.icvs)) 1:4 else 1:3
    legend("topright", title=analog, 
           legend=c("Analog", "Target", "Separation of means", "ICV")[s], 
           pt.bg=c("dodgerblue", "grey", NA, NA)[s], 
           col = c("black", "black", "black", "black")[s], 
           pt.cex=c(1.6,1.3,1.3, 1)[s], 
           pch=c(21, 21, 17, 16)[s], 
           bty="n")
    box()
  }
  
  ## 2D scatterplot
  if(plot2d){
    par(mfrow=c(2,2), mar=c(3,3,1,1), mgp=c(1.75,0.25,0))
    for(i in 1:4){
      a <- predict(pca, clim.analog)[, plot2d.pcs[i,]]
      b <- predict(pca, clim.target)[, plot2d.pcs[i,]]
      b <- sweep(b, 2, apply(a, 2, mean), '-') # shift the target data so that the analog centroid is at zero. this is done at a later stage than the pca in the distance calculation.
      a <- sweep(a, 2, apply(a, 2, mean), '-') # centre the analog centroid on zero. this is done at a later stage than the pca in the distance calculation.
      plot(a, col="dodgerblue", xlim=range(c(a[,1], b[,1])), ylim=range(c(a[,2], b[,2])), asp=1, tck=0.01)
      points(b, bg=ColScheme[cut(q, breakpoints)], pch=21, cex=1.5)
      if(!is.null(clim.icvs)){
        c <- predict(pca, clim.icv)[, plot2d.pcs[i,]]
        c <- sweep(c, 2, apply(c, 2, mean), '-') # centre the ICV on the analog centroid. this is done at a later stage than the pca in the distance calculation. 
        points(c, col="black", pch=16, cex=1)
      }
      points(a, col="dodgerblue", pch=16)
      mtext(paste(analog, "\n", pcs, "PCs"), line=-2.5, adj = 0.05, )
    }
  }
  
  ## 3D scatterplot
  if(plot3d){
    
    # revert to the raw pcs (centered on the analog centroid), because standardization obscures the shape of the analog distribution
    a <- predict(pca, clim.analog)
    b <- predict(pca, clim.target)
    b <- sweep(b, 2, apply(a, 2, mean), '-') # shift the target data so that the analog centroid is at zero. this is done at a later stage than the pca in the distance calculation.
    a <- sweep(a, 2, apply(a, 2, mean), '-') # centre the analog centroid on zero. this is done at a later stage than the pca in the distance calculation.
    
    b_colors <- ColScheme[cut(q, breakpoints)] # Define colors for points in 'b'
    
    # Create the 3D scatterplot
    plot <- plot_ly() %>%
      add_trace(
        x = a[, plot3d.pcs[1]], y = a[, plot3d.pcs[2]], z = a[, plot3d.pcs[3]],
        type = "scatter3d", mode = "markers",
        marker = list(size = 5, color = "dodgerblue", opacity = 1),
        name = "Analog Points"
      ) %>%
      add_trace(
        x = b[, plot3d.pcs[1]], y = b[, plot3d.pcs[2]], z = b[, plot3d.pcs[3]],
        type = "scatter3d", mode = "markers",
        marker = list(size = 6, color = b_colors, opacity = 1),
        name = "Target Points"
      ) 
    # Add ICV points if they exist
    if(!is.null(clim.icvs)) {
      c <- predict(pca, clim.icv)
      c <- sweep(c, 2, apply(c, 2, mean), '-') # centre the ICV on the analog centroid. this is done at a later stage than the pca in the distance calculation. 
      plot <- plot %>%
        add_trace(
          x = c[, plot3d.pcs[1]], y = c[, plot3d.pcs[2]], z = c[, plot3d.pcs[3]],
          type = "scatter3d", mode = "markers",
          marker = list(size = 4, color = "black", opacity = 1),
          name = "ICV"
        )
    }
    # Add biplot lines
    if(biplot) {
      loadings <- pca$rotation[, plot3d.pcs]
      scale_factor <- max(abs(c(a, b))) * 2
      scaled_loadings <- loadings * scale_factor
      for (i in 1:nrow(scaled_loadings)) {
        plot <- plot %>%
          add_trace(
            x = c(0, scaled_loadings[i, 1]),
            y = c(0, scaled_loadings[i, 2]),
            z = c(0, scaled_loadings[i, 3]),
            type = "scatter3d",
            mode = "lines+text",
            line = list(color = "black", width = 2),
            text = rownames(scaled_loadings)[i],
            textposition = "middle center",
            showlegend = FALSE, 
            name = paste("Loading:", rownames(scaled_loadings)[i])
          )
      }
    }
    plot <- plot %>%
      layout(
        scene = list(
          xaxis = list(title = paste0("PC", plot3d.pcs[1])),
          yaxis = list(title = paste0("PC", plot3d.pcs[2])),
          zaxis = list(title = paste0("PC", plot3d.pcs[3]))
        ),
        title = list(text = paste(analog, "\nNovelty in", pcs, "PCs"), x = 0.05)
      )
    # Display the plot
    print(plot)
  }
  return(novelty)
}

# gcm_weight <- data.table(gcm = c("ACCESS-ESM1-5", "BCC-CSM2-MR", "CanESM5", "CNRM-ESM2-1", "EC-Earth3", 
#                                  "GFDL-ESM4", "GISS-E2-1-G", "INM-CM5-0", "IPSL-CM6A-LR", "MIROC6", 
#                                  "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL"),
#                          weight = c(1,0,0,1,1,1,1,0,0,1,1,1,0))
# 
# rcp_weight <- data.table(rcp = c("ssp126","ssp245","ssp370","ssp585"), 
#                          weight = c(0.8,1,0.8,0))
# 
# all_weight <- as.data.table(expand.grid(gcm = gcm_weight$gcm,rcp = rcp_weight$rcp))
# all_weight[gcm_weight,wgcm := i.weight, on = "gcm"]
# all_weight[rcp_weight,wrcp := i.weight, on = "rcp"]
# all_weight[,weight := wgcm*wrcp]
# all_weight[,comb := paste0("('",gcm,"','",rcp,"',",weight,")")]
# weights <- paste(all_weight$comb,collapse = ",")
# 
# temp <- dbGetQuery(conn, q2)
# setDT(temp)
# t2 <- unique(temp[,.(siteno,gcm,scenario,run)])
# 
# siteno <- c(43,47,48,49,51)
# 
# groupby = "siteno"
# cciss_sql <- paste0("
#   WITH cciss AS (
#     SELECT cciss_future13_array.siteno,
#          labels.gcm,
#          labels.scenario,
#          labels.futureperiod,
#          labels.run,
#          bgc_attribution13.bgc,
#          bgcv13.bgc bgc_pred,
#          w.weight
#   FROM cciss_future13_array
#   JOIN bgc_attribution13
#     ON (cciss_future13_array.siteno = bgc_attribution13.siteno),
#        unnest(bgc_id) WITH ordinality as source(bgc_id, row_idx)
#   JOIN (SELECT ROW_NUMBER() OVER(ORDER BY gcm_id, scenario_id, futureperiod_id, run_id) row_idx,
#                gcm,
#                scenario,
#                futureperiod,
#                run
#         FROM gcm 
#         CROSS JOIN scenario
#         CROSS JOIN futureperiod
#         CROSS JOIN run) labels
#     ON labels.row_idx = source.row_idx
#     JOIN (values ",weights,") 
#     AS w(gcm,scenario,weight)
#     ON labels.gcm = w.gcm AND labels.scenario = w.scenario
#   JOIN bgcv13
#     ON bgcv13.bgc_id = source.bgc_id
#   WHERE cciss_future13_array.siteno IN (", paste(unique(siteno), collapse = ","), ")
#   AND futureperiod IN ('2001', '2021','2041','2061','2081')
#   
#   ), cciss_count_den AS (
#   
#     SELECT ", groupby, " siteref,
#            futureperiod,
#            SUM(weight) w
#     FROM cciss
#     GROUP BY ", groupby, ", futureperiod
#   
#   ), cciss_count_num AS (
#   
#     SELECT ", groupby, " siteref,
#            futureperiod,
#            bgc,
#            bgc_pred,
#            SUM(weight) w
#     FROM cciss
#     GROUP BY ", groupby, ", futureperiod, bgc, bgc_pred
#   
#   ), cciss_curr AS (
#       SELECT cciss_prob13.siteno,
#       '1991' as period,
#       bgc_attribution13.bgc,
#       bgc_pred,
#       prob
#       FROM cciss_prob13
#       JOIN bgc_attribution13
#       ON (cciss_prob13.siteno = bgc_attribution13.siteno)
#       WHERE cciss_prob13.siteno IN (", paste(unique(siteno), collapse = ","), ")
#       
#   ), curr_temp AS (
#     SELECT ", groupby, " siteref,
#            COUNT(distinct siteno) n
#     FROM cciss_curr
#     GROUP BY ", groupby, "
#   )
#   
#   SELECT cast(a.siteref as text) siteref,
#          a.futureperiod,
#          a.bgc,
#          a.bgc_pred,
#          a.w/cast(b.w as float) bgc_prop
#   FROM cciss_count_num a
#   JOIN cciss_count_den b
#     ON a.siteref = b.siteref
#    AND a.futureperiod = b.futureperiod
#    WHERE a.w <> 0
#   
#   UNION ALL
# 
#   SELECT cast(", groupby, " as text) siteref,
#           period as futureperiod,
#           bgc,
#           bgc_pred,
#           SUM(prob)/b.n bgc_prop
#   FROM cciss_curr a
#   JOIN curr_temp b
#     ON a.",groupby," = b.siteref
#   WHERE siteno in (", paste(unique(siteno), collapse = ","), ")
#   GROUP BY ", groupby, ",period,b.n, bgc, bgc_pred
#   
#   UNION ALL
# 
#   SELECT DISTINCT 
#             cast(", groupby, " as text) siteref,
#             '1961' as futureperiod,
#             bgc,
#             bgc as bgc_pred,
#             cast(1 as numeric) bgc_prop
#     FROM cciss_curr
#     WHERE siteno IN (", paste(unique(siteno), collapse = ","), ")
#   ")
# 
# 
# test <- dbGetQuery(conn, cciss_sql)
# setDT(test)
# setorder(test, siteref, futureperiod)
# 
