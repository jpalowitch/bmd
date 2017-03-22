library(RColorBrewer)
library(plotrix)
library(igraph)

total_expers <- readLines("sims-results/exper-names.txt")

source("sim_eQTL_network.R")
source("best_match.R")

# Set method names:
methNames <- c("bmd", "bmd2")

# Set which methods to plot and their plot names
plot_meths <- c(1:2)
plot_names <- c("BMD", "BMD2")

# Set points
pchs <- c(14, 8, 3, 22, 24, 21)

# Re-get results?
getResults <- TRUE

# Plot main text?
main_text_plot <- TRUE

plot_expers <- 9

# This should consistent throughout the experiments
nreps <- 20

# Brewing colors
colPal <- brewer.pal(9, "Set1")


for (exper in plot_expers) {
  
  expString <- paste0("experiment", total_expers[exper])
  
  # Finding the folder
  root_dir <- file.path("sims-results", expString)
  
  # Getting par values
  if (exists("axis_par_string")) {rm("axis_par_string")}
  load(file.path("sims-results/sbm-par-lists", paste0(expString, ".RData")))
  
  allscores <- array(0, dim = c(par_divs, nreps, 8))
  methscores <- rep(list(allscores), length(methNames))
  names(methscores) <- methNames
  
  if (getResults) {
    
    for (p in 1:par_divs) {
      
      cat("########\n")
      cat("p=",p,"\n")
      cat("########\n")
      
      curr_dir_p <- file.path(root_dir, par_dirs[p])
      
      for (rep in 1:nreps) {
          
        cat("rep =",rep,"\n")
        curr_dir_p_rep <- file.path(curr_dir_p, rep)
        load(file.path(curr_dir_p_rep, "sim.RData"))
        dx <- ncol(sim$X); dy <- ncol(sim$Y); n <- nrow(sim$X)
                
        for (meth in methNames) {
            
          meth_fn <- file.path(curr_dir_p_rep, paste0(meth, ".RData"))
          load(meth_fn)
          
          comms <- results$communities
          full_comms <- lapply(seq_along(comms$X_sets),
                               function (i) c(comms$X_sets[[i]],                                             comms$Y_sets[[i]]))
          bg <- setdiff(1:(dx + dy), 
                        c(unlist(comms$X_sets), unlist(comms$Y_sets)))
          true_bg <- setdiff(1:(dx + dy),
                             c(unlist(sim$bm)))
          scores <- best_match_bimodule(full_comms, sim$bms,
                                        bg, true_bg)
          methscores[[meth]][p, rep, ] <- scores
        
          rm(results)
          gc()
          
        }
        
        rm(sim)
        gc()
        
      }
      
    }
    
    
      
    convertToMat <- function (resultsList, type = "mean", score = 1) {
      
      if (type == "mean") {
        stats <- unlist(lapply(resultsList, 
                               function(arr) apply(arr[ , , score], 1, mean)))
      } else {
        stats <- unlist(lapply(resultsList, 
                               function(arr) apply(arr[ , , score], 1, sd)))
      }
      
      stats <- matrix(stats, nrow = length(methNames), byrow = T)
      rownames(stats) <- methNames
      return(stats)
      
    }
      
    BM_means <- convertToMat(methscores, score = 1)
    BM1_means <- convertToMat(methscores, score = 2)
    BM2_means <- convertToMat(methscores, score = 3)
    BJ_means <- convertToMat(methscores, score = 4)
    SP_means <- convertToMat(methscores, score = 5)
    SM_means <- convertToMat(methscores, score = 6)
    F1_means <- convertToMat(methscores, score = 7)
    F2_means <- convertToMat(methscores, score = 8)
      
    BM_sds <- convertToMat(methscores, score = 1, type = "sd")
    BM1_sds <- convertToMat(methscores, score = 2, type = "sd")
    BM2_sds <- convertToMat(methscores, score = 3, type = "sd")
    BJ_sds <- convertToMat(methscores, score = 4, type = "sd")
      
    save(BM_means, BM_sds,
         BM1_means, BM1_sds,
         BM2_means, BM2_sds,
         BJ_means, BJ_sds,
         SP_means,
         SM_means,
         F1_means,
         F2_means,
         file = file.path(root_dir, "plot_results.RData"))
    
  } else {
    
    load(file.path(root_dir, "plot_results.RData"))
    
  }
  
  # Color palette
  colPal <- colPal[1:length(methNames)]
  
  source("makePerformancePlot.R")
  
  # Some plot defaults
  cex.main <- 2.5
  cex.lab <- 2.5
  cex.axis <- 2.5
  cex <- 3.5
  legCex <- 2.5
  lwd <- 3
  dotnmi <- FALSE
  
  if (exists("axis_par_string")) {
    xlab_string <- axis_par_string
  } else {
    xlab_string <- pars[axis_par]
  }
  
  paramVec <- par_settings[axis_par, ]
  
  if (paste0(expString, "_BM.png") %in% 
      list.files(file.path("sims-results", expString))) {
    file.remove(paste0(expString, "_BM.png"))
  }

  if (main_text_plot) {
    main_str <- main_text
  } else {
    main_str <- ""
  }
  
  
  BM_means <- BM_means[plot_meths, , drop = FALSE]
  BM_sds <- BM_sds[plot_meths, , drop = FALSE]
  rownames(BM_means) <- plot_names
  rownames(BM_sds) <- plot_names
  
  BM1_means <- BM1_means[plot_meths, , drop = FALSE]
  BM1_sds <- BM1_sds[plot_meths, , drop = FALSE]
  rownames(BM1_means) <- plot_names
  rownames(BM1_sds) <- plot_names
  
  BM2_means <- BM2_means[plot_meths, , drop = FALSE]
  BM2_sds <- BM2_sds[plot_meths, , drop = FALSE]
  rownames(BM2_means) <- plot_names
  rownames(BM2_sds) <- plot_names
  
  BJ_means <- BJ_means[plot_meths, , drop = FALSE]
  BJ_sds <- BJ_sds[plot_meths, , drop = FALSE]
  rownames(BJ_means) <- plot_names
  rownames(BJ_sds) <- plot_names
  
  SP_means <- SP_means[plot_meths, , drop = FALSE]
  SM_means <- SM_means[plot_meths, , drop = FALSE]
  F1_means <- F1_means[plot_meths, , drop = FALSE]
  F2_means <- F2_means[plot_meths, , drop = FALSE]
  rownames(SP_means) <- rownames(SM_means) <- rownames(F1_means) <- 
    rownames(F2_means) <- plot_names
  
  
  plotfn = file.path("sims-results", paste0(expString, ".png"))
  
  png(plotfn, width = 1500, height = 1500)
  par(mfrow = c(3, 3), oma = rep(0, 4),
      mar = c(11, 11, 6, 4),
      mgp = c(6, 2, 0))
  
  # BM
  
  suppressWarnings(
    
    dummy <- makePerformancePlot(plotFile = FALSE, 
                                 xvals = paramVec, yRange = c(0, 1),
                                 meanMat = BM_means, tnmi = dotnmi,
                                 sdMat = BM_sds,
                                 xRange = c(paramVec[1], paramVec[length(paramVec)]),
                                 main = main_str,
                                 xlab = xlab_string,
                                 ylab = "Best Match Jaccard",
                                 legPos = "bottomright",
                                 legCex = legCex,
                                 lwd = lwd,
                                 cex = cex, pchs = pchs,
                                 cex.main = cex.main,
                                 cex.lab = cex.lab,
                                 cex.axis = cex.axis)
  )
    
  # BM1
  
  suppressWarnings(
    
    dummy <- makePerformancePlot(plotFile = FALSE, doLegend = FALSE,
                                 xvals = paramVec, yRange = c(0, 1),
                                 meanMat = BM1_means, tnmi = dotnmi,
                                 sdMat = BM1_sds,
                                 xRange = c(paramVec[1], paramVec[length(paramVec)]),
                                 main = main_str,
                                 xlab = xlab_string,
                                 ylab = "Best Match Res.prop",
                                 legPos = "bottomright",
                                 legCex = legCex,
                                 lwd = lwd,
                                 cex = cex, pchs = pchs,
                                 cex.main = cex.main,
                                 cex.lab = cex.lab,
                                 cex.axis = cex.axis)
    
  )
  
  # BM2

  suppressWarnings(
    
    dummy <- makePerformancePlot(plotFile = FALSE, doLegend = FALSE,
                                 xvals = paramVec, yRange = c(0, 1),
                                 meanMat = BM2_means, tnmi = dotnmi,
                                 sdMat = BM2_sds,
                                 xRange = c(paramVec[1], paramVec[length(paramVec)]),
                                 main = main_str,
                                 xlab = xlab_string,
                                 ylab = "Best Match Truth.prop",
                                 legPos = "topright",
                                 legCex = legCex,
                                 lwd = lwd,
                                 cex = cex, pchs = pchs,
                                 cex.main = cex.main,
                                 cex.lab = cex.lab,
                                 cex.axis = cex.axis)
    
  )
  
  # BJ

  suppressWarnings(
    
    dummy <- makePerformancePlot(plotFile = FALSE, doLegend = FALSE,
                                 xvals = paramVec, yRange = c(0, 1),
                                 meanMat = BJ_means, tnmi = dotnmi,
                                 sdMat = BJ_sds,
                                 xRange = c(paramVec[1], paramVec[length(paramVec)]),
                                 main = main_str,
                                 xlab = xlab_string,
                                 ylab = "Background Jaccard",
                                 legPos = "topright",
                                 legCex = legCex,
                                 lwd = lwd,
                                 cex = cex, pchs = pchs,
                                 cex.main = cex.main,
                                 cex.lab = cex.lab,
                                 cex.axis = cex.axis)
    
  )
  
  # StickyProb
  
  suppressWarnings(
    
    dummy <- makePerformancePlot(plotFile = FALSE, doLegend = FALSE,
                                 xvals = paramVec, yRange = c(0, 1),
                                 meanMat = SP_means, tnmi = dotnmi,
                                 sdMat = NA,
                                 xRange = c(paramVec[1], paramVec[length(paramVec)]),
                                 main = main_str,
                                 xlab = xlab_string,
                                 ylab = "StickyProb",
                                 legPos = "topright",
                                 legCex = legCex,
                                 lwd = lwd,
                                 cex = cex, pchs = pchs,
                                 cex.main = cex.main,
                                 cex.lab = cex.lab,
                                 cex.axis = cex.axis)
    
  )
  
  # StickyMean
  
  suppressWarnings(
    
    dummy <- makePerformancePlot(plotFile = FALSE, doLegend = FALSE,
                                 xvals = paramVec, yRange = c(0, 1),
                                 meanMat = SP_means, tnmi = dotnmi,
                                 sdMat = NA,
                                 xRange = c(paramVec[1], paramVec[length(paramVec)]),
                                 main = main_str,
                                 xlab = xlab_string,
                                 ylab = "StickyMean",
                                 legPos = "topright",
                                 legCex = legCex,
                                 lwd = lwd,
                                 cex = cex, pchs = pchs,
                                 cex.main = cex.main,
                                 cex.lab = cex.lab,
                                 cex.axis = cex.axis)
    
  )
  
  # Avg FDR
  
  suppressWarnings(
    
    dummy <- makePerformancePlot(plotFile = FALSE, doLegend = FALSE,
                                 xvals = paramVec, yRange = c(0, 1),
                                 meanMat = SP_means, tnmi = dotnmi,
                                 sdMat = NA,
                                 xRange = c(paramVec[1], paramVec[length(paramVec)]),
                                 main = main_str,
                                 xlab = xlab_string,
                                 ylab = "Avg FDR",
                                 legPos = "topright",
                                 legCex = legCex,
                                 lwd = lwd,
                                 cex = cex, pchs = pchs,
                                 cex.main = cex.main,
                                 cex.lab = cex.lab,
                                 cex.axis = cex.axis)
    
  )
  
  # AvgWtdFDR
  
  suppressWarnings(
    
    dummy <- makePerformancePlot(plotFile = FALSE, doLegend = FALSE,
                                 xvals = paramVec, yRange = c(0, 1),
                                 meanMat = SP_means, tnmi = dotnmi,
                                 sdMat = NA,
                                 xRange = c(paramVec[1], paramVec[length(paramVec)]),
                                 main = main_str,
                                 xlab = xlab_string,
                                 ylab = "AvgWtdFDR",
                                 legPos = "topright",
                                 legCex = legCex,
                                 lwd = lwd,
                                 cex = cex, pchs = pchs,
                                 cex.main = cex.main,
                                 cex.lab = cex.lab,
                                 cex.axis = cex.axis)
    
  )
  
  dev.off()
  
}