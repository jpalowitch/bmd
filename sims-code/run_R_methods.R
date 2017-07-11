#run_expers <- sapply(commandArgs(TRUE), as.numeric)
run_expers <- 15:17
source("sim_eQTL_network.R")
source("run_brim.R")
source("ircc.R")
source("cbceNW.R")
source("cbceNW_c.R")
total_expers <- readLines("sims-results/exper-names.txt")

runBMDcpp <- TRUE
runBMD2 <- TRUE
runBMD_c <- TRUE
runBMD <- FALSE
runBRIM <- FALSE
runkmeans <- FALSE

# This should consistent throughout the experiments
# (and match the same variable in sims/lfr/make_lfr_sims.R)
nreps <- 10

for (exper in run_expers) {
  
  exper_string <- paste0("experiment", total_expers[exper])
  
  # Finding the folder
  root_dir <- file.path("sims-results", exper_string)
  
  # Loading parameters
  load(paste0(file.path("sims-results/sbm-par-lists", exper_string), ".RData"))
  
  for (p in 1:par_divs) {
    
    curr_dir_p <- file.path(root_dir, par_dirs[p])
    
    # Setting alpha
    alpha <- ifelse(palpha, par_settings[1, p], 0.05)
    
    for (rep in 1:nreps) {
      
      cat("exper", exper, "p", p, "rep", rep, "\n")
    
      curr_dir_p_rep <- file.path(curr_dir_p, rep)
      load(file.path(curr_dir_p_rep, "sim.RData"))
      
      # Running BMDcpp
      if (runBMD_c) {
        timer <- proc.time()[3]
        results <- cbceNW_c(sim$X, sim$Y, alpha = alpha, verbose = TRUE, generalOutput = TRUE,
                            updateOutput = TRUE, OL_tol = 100, Dud_tol = Inf, Cpp = FALSE, pval_parallel = TRUE,
                            calc_full_cor = TRUE, updateMethod = 1, twoSided = TRUE, parallel = FALSE)
        timer <- proc.time()[3] - timer
        save(results, timer, file = file.path(curr_dir_p_rep, "bmd_c.RData"))
      }
      
      # Running BMDcpp
      if (runBMDcpp) {
        timer <- proc.time()[3]
        results <- cbceNW(sim$X, sim$Y, alpha = alpha, verbose = FALSE, generalOutput = TRUE,
                          updateOutput = FALSE, OL_tol = 100, Dud_tol = Inf, Cpp = TRUE,
                          calc_full_cor = TRUE, updateMethod = 1)
        timer <- proc.time()[3] - timer
        save(results, timer, file = file.path(curr_dir_p_rep, "bmd_cpp.RData"))
      }
      
      # Running BMD2
      if (runBMD2) {
        timer <- proc.time()[3]
        results <- cbceNW(sim$X, sim$Y, alpha = alpha, verbose = FALSE, generalOutput = TRUE,
                          updateOutput = FALSE, OL_tol = 100, Dud_tol = Inf,
                          calc_full_cor = TRUE, updateMethod = 2)
        timer <- proc.time()[3] - timer
        save(results, timer, file = file.path(curr_dir_p_rep, "bmd2.RData"))
      }
      
      # Running BMD
      if (runBMD) {
        timer <- proc.time()[3]
        results <- cbceNW(sim$X, sim$Y, alpha = alpha, verbose = FALSE,
                          updateOutput = FALSE, OL_tol = 100, Dud_tol = 50,
                          calc_full_cor = TRUE, updateMethod = 1)
        timer <- proc.time()[3] - timer
        save(results, timer, file = file.path(curr_dir_p_rep, "bmd.RData"))
      }
      
      # Running BRIM
      if (runBRIM) {
        timer <- proc.time()[3]
        results <- run_brim(sim$X, sim$Y, alpha = alpha)
        timer <- proc.time()[3] - timer
        save(results, timer, file = file.path(curr_dir_p_rep, "brim.RData"))
      }
      
      # Running IRCC-kmeans
      if (runkmeans) {
        nbmds <- ifelse(par_list$bgmult > 0, par_list$b + 1, par_list$b)
        timer <- proc.time()[3]
        results <- ircc(sim$X, sim$Y, nbmds = nbmds, method = "kmeans")
        timer <- proc.time()[3] - timer
        save(results, timer, file = file.path(curr_dir_p_rep, "kmeans.RData"))
      }
        
      rm(sim, results)
      gc()
      
    }
    
  }
  
}