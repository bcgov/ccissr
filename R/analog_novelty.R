
# Copyright 2024 Province of British Columbia
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' Novelty (outlyingness) measurement for pre-determined spatial climate analogs. 
#'
#' @description
#' 
#' This function calculates novelty (or outlyingness) of climate conditions by measuring 
#' the Mahalanobis distance between these target climates and the centroid of their pre-determined
#' climate analogs. The covariance matrix of the Mahalanobis distance is the spatial 
#' and/or interannual variation in climate conditions of the climate analog, as 
#' selected by the user. It optionally renders 2D and 3D plots for visualizing 
#' the distance measurements in the principal component space.
#' 
#' @details
#' 
#' The novelty measure is calculated using Mahalanobis distance, which quantifies the 
#' statistical distance between the target climate conditions and their spatial analogs 
#' (i.e., climates similar to the target in historical records). The function can also 
#' incorporate interannual climatic variability (ICV) in the covariance matrix of the climate 
#' analog to provide additional information on the climatic scale of each principal component. 
#' For visualization purposes, the function supports 2D and 3D scatterplots 
#' based on principal component analysis (PCA), showing how the target climate conditions 
#' relate to the analogs and, optionally, to the ICV.
#' 
#' The main output of the function is a vector of novelty values, either as Mahalanobis 
#' distances or as sigma dissimilarities (Interpretable as the number of standard deviations away from 
#' the analog cluster). These values are returned for each target climate in the dataset. 
#' Additional outputs include several optional plots:
#' - A scree plot (if `plotScree` is enabled) showing the standard deviation of the 
#' analog and target variation in the principal components.
#' - A 2D scatterplot (if `plot2d` is enabled) showing the relationship between the target 
#'   and analogs across selected principal components
#' - A 3D scatterplot (if `plot3d` is enabled) to visualize the distribution of the 
#'   target and analogs in 3D PCA space, with optional biplot vectors representing the 
#'   correlation of the standardized climate variables with the principal components.
#' 
#' The function supports customization of the principal component selection (via `threshold` or 
#' `pcs`), weighting of the ICV in the distance measurement (`weight.icv`), and the selection 
#' of specific analogs for visualization (`analog.focal`).
#' 
#'  
#' @param clim.targets data.table. Climate variables for which the analogs were 
#' identified
#' @param clim.analogs data.table. Climate variables of at least 50 locations 
#' representing the spatial variation in the reference period (historical) climate
#' of each analog in the analog pool. 
#' @param label.targets character. Vector of the analog IDs identified for the 
#' climatic conditions listed in `clim.targets`. Length equals number of records 
#' in `clim.targets`.  
#' @param label.analogs character. Vector of the analog IDs for the climatic 
#' conditions listed in `clim.analogs`. Length equals number of records in
#' `clim.analogs`.  
#' @param vars character. Climate variables to use in the novelty measurement. 
#' @param clim.icvs data.table. Time series of climate variables at the geographic
#' centroids of each analog in the analog pool. If not null, this interannual climatic
#' variability will be pooled with the spatial variation of the analog to calculate 
#' the covariance matrix used in the Mahalanobis distance measurement. 
#' @param label.icvs character. Vector of the analog IDs for the climatic 
#' conditions listed in `clim.icvs`. Length equals number of records in `clim.icvs`.  
#' @param analog.focal character. Optionally specify a single analog for visualization. 
#' @param threshold numeric. The cumulative variance explained to use as a threshold
#' for truncation of principal components to use in the mahalanobis distance measurement. 
#' @param pcs integer. The fixed number of PCs to use in the mahalanobis distance 
#' measurement. Non-null values of this parameter override the `threshold` parameter.
#' @param logVars logical. Log-transform ratio variables, i.e., variables with a minimum 
#' of zero. 
#' @param plotScree logical. If `analog.focal` is specified, plot a scree plot showing 
#' the standard deviation of the analog points, target points, and (if activated) interannual
#' climatic variability (ICV). The difference of the analog and target means is also shown. 
#' @param plot2d logical. If `analog.focal` is specified, plot a 4-panel set of 
#' bivariate scatterplots visualizing the distances in the principal components. 
#' @param plot2d.pcs numeric matrix. a 4x2 matrix indicating the principal components to 
#' display in the 2D scatterplots. Rows correspond to the panel to be displayed, 
#' and columns indicate the principal components to display in the x and y axes. 
#' @param plot3d logical. If `analog.focal` is specified, plot a 3-dimensional scatterplot 
#' visualizing the distances in the first three PCs. 
#' @param plot3d.pcs numeric. Principal components to display in the x, y, and z 
#' axes of the 3D scatterplot. 
#' @param biplot logical. Include lines on the 3D plot indicating the correlation 
#' of the standardized climate variables with the principal components. 
#' @param plot3d.candidates logical. Include representative points and labelled centroids for all 
#' western North American candidate analogs. 
#' 
#' @return `vector` of sigma dissimilarity (if sigma==TRUE) or Mahalanobis distances 
#' (if sigma==FALSE) corresponding to each element of the 'analogs.target' vector. 
#'
#' @importFrom stats prcomp mahalanobis cov cor sd
#' @importFrom graphics plot text points abline arrows
#' @importFrom grDevices rgb
#' @importFrom plotly plot_ly layout add_trace
#' @importFrom EnvStats pchi qchi
#'  
#' @examples
#' if (FALSE) {
#'   
#' }
#' #'
#' @export

analog_novelty <- function(clim.targets, clim.analogs, label.targets, label.analogs, vars,
                           clim.icvs = NULL, label.icvs = NULL, weight.icv = 0.5, sigma = TRUE,
                           analog.focal = NULL, threshold = 0.95, pcs = NULL, logVars = TRUE, 
                           plotScree = FALSE, 
                           plot2d = FALSE, plot2d.pcs = cbind(c(1,2,3,4), c(2,3,4,5)), 
                           plot3d = FALSE, plot3d.pcs=c(1,2,3), biplot = TRUE, 
                           plot3d.candidates = FALSE){
  
  analogs <- if(is.null(analog.focal)) unique(label.targets) else analog.focal # list of analogs to loop through
  novelty <- rep(NA, length(label.targets)) # initiate a vector to store the sigma dissimilarities
  
  for(analog in analogs){ # loop through all of the analogs used to describe the target climates. 
    clim.analog <- clim.analogs[label.analogs==analog, ..vars]
    clim.target <- clim.targets[label.targets==analog, ..vars]
    if(!is.null(clim.icvs)) clim.icv <- clim.icvs[label.icvs==analog, ..vars]
    if(plot3d.candidates) clim.analogs.all <- clim.analogs[, ..vars]
      
    ## data cleaning
    clim.analog <- clim.analog[complete.cases(clim.analog)] # remove rows without data
    clim.analog <- clim.analog[, .SD, .SDcols = which(sapply(clim.analog, function(x) var(x, na.rm = TRUE) > 0))]  # Remove zero-variance columns
    clim.target <- clim.target[, .SD, .SDcols = names(clim.analog)]
    if(!is.null(clim.icvs)) clim.icv <- clim.icv[complete.cases(clim.icv)]
    if(!is.null(clim.icvs)) clim.icv <- clim.icv[, .SD, .SDcols = names(clim.analog)]
    if(plot3d.candidates){
      clim.analogs.all <- clim.analogs.all[complete.cases(clim.analogs.all)]
    clim.analogs.all <- clim.analogs.all[, .SD, .SDcols = names(clim.analog)]
    }
    
    ## log-transform ratio variables
    if(logVars){
      clim.analog <- logVars(clim.analog, zero_adjust = TRUE)
      clim.target <- logVars(clim.target, zero_adjust = TRUE)
      if(!is.null(clim.icvs)) clim.icv <- logVars(clim.icv, zero_adjust = TRUE)
      if(plot3d.candidates) clim.analogs.all <- logVars(clim.analogs.all, zero_adjust = TRUE)
      
      ## remove variables with non-finite values in the target population (this is an edge case that occurs when the target population has a variable (typically CMD) with only zeroes)
      clim.target <- clim.target[, lapply(.SD, function(x) if (all(is.finite(x))) x else NULL)]
      clim.analog <- clim.analog[, .SD, .SDcols = names(clim.target)]
      if(!is.null(clim.icvs)) clim.icv <- clim.icv[, .SD, .SDcols = names(clim.target)]
      if(plot3d.candidates) clim.analogs.all <- clim.analogs.all[, .SD, .SDcols = names(clim.target)]
    }
    
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
    if(plot3d.candidates) clim.analogs.all[, (names(clim.analogs.all)) := lapply(names(clim.analogs.all), function(col) {
      (get(col) - unlist(clim.mean)[col]) / unlist(clim.sd)[col]
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
    if(plot3d.candidates){
      d <- predict(pca, clim.analogs.all)
      d <- sweep(d, 2, apply(a, 2, mean), '-') # shift the candidate data so that the analog centroid is at zero. 
      e <- aggregate(d, by=list(label.analogs), FUN=mean)
      e <- e[is.finite(e$PC1),]
      label.analogs.mean <- e$Group.1
      e <- e[,-1]
    }
    a <- sweep(a, 2, apply(a, 2, mean), '-') # centre the analog centroid on zero. this is done at a later stage than the pca in the distance calculation.
    
    b_colors <- ColScheme[cut(q, breakpoints)] # Define colors for points in 'b'
    
    # Create the 3D scatterplot
    plot <- plot_ly() %>%
      add_trace(
        x = a[, plot3d.pcs[1]], y = a[, plot3d.pcs[2]], z = a[, plot3d.pcs[3]],
        type = "scatter3d", mode = "markers",
        marker = list(size = 5, color = "dodgerblue", opacity = 1),
        hoverinfo = "none", # Turn off hover labels
        name = "Analog Points"
      ) %>%
      add_trace(
        x = b[, plot3d.pcs[1]], y = b[, plot3d.pcs[2]], z = b[, plot3d.pcs[3]],
        type = "scatter3d", mode = "markers",
        marker = list(size = 6, color = b_colors, opacity = 1),
        hoverinfo = "none", # Turn off hover labels
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
          hoverinfo = "none", # Turn off hover labels
          name = "ICV"
        )
    }
    # Add candidate analogs
    if(plot3d.candidates){
      data("zones_colours_ref")
      zone <- rep(NA, length(label.analogs.mean))
      for(i in zones_colours_ref$classification){ zone[grep(i,label.analogs.mean)] <- i }
      # zone <- factor(zone, zones_colours_ref$classification)
      zone_colours <- as.character(zones_colours_ref$colour[match(zone, zones_colours_ref$classification)]) 
      zone_colours[is.na(zone_colours)] <- "#808080"
      # zone_colours <- factor(zone_colours, zones_colours_ref$colour)
      
      plot <- plot %>%
        add_trace(
          x = d[, plot3d.pcs[1]], y = d[, plot3d.pcs[2]], z = d[, plot3d.pcs[3]],
          type = "scatter3d", mode = "markers",
          marker = list(size = 2, color = "#cccccc", opacity = 0.35),
          hoverinfo = "none", # Turn off hover labels
          name = "All analogs"
        ) %>%
        add_trace(
          x = e[, plot3d.pcs[1]], y = e[, plot3d.pcs[2]], z = e[, plot3d.pcs[3]],
          type = "scatter3d", mode = "markers+text",
          marker = list(size = 5, color = zone_colours, opacity = 1),
          # marker = list(size = 3, color = "#000000", opacity = 0.5),
          text = label.analogs.mean, # Vector of labels corresponding to points in e
          textposition = "right",
          textfont = list(size = 8, color = "#666666"),
          hoverinfo = "none", # Turn off hover labels
          name = "BGC centroids"
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
          xaxis = list(title = paste0("PC", plot3d.pcs[1]), showspikes = FALSE),
          yaxis = list(title = paste0("PC", plot3d.pcs[2]), showspikes = FALSE),
          zaxis = list(title = paste0("PC", plot3d.pcs[3]), showspikes = FALSE)
        ),
        title = list(text = paste(analog, "\nNovelty in", pcs, "PCs"), x = 0.05)
      )
    # Display the plot
    print(plot)
  }
  return(novelty)
}

