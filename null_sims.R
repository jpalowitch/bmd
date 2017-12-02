source("mvrnormR.R")
library(cbce)

nsims <- 20
ns <- seq(50, 500, 50)
rhos <- seq(0, 0.9, 0.1)
dx <- 100
mx <- 50
dy <- 100
my <- 50

# ------------------------------------------------------------------------------

# Setting up 
SigmaX <- matrix(0, dx, dx)
SigmaY <- matrix(0, dy, dy)
meth1_ncomm <- meth2_ncomm <- diffs_ncomm <- matrix(0, length(ns), length(rhos))
meth1_nnode <- meth2_nnode <- diffs_nnode <- matrix(0, length(ns), length(rhos))
I <- length(ns)
J <- length(rhos)
K <- nsims

# ------------------------------------------------------------------------------

for (i in seq_along(ns)) {
  
  n <- ns[i]
  
  for (j in seq_along(rhos)) {
    
    rho <- rhos[j]
    SigmaX[1:mx, 1:mx] <- SigmaY[1:mx, 1:mx] <- rho
    diag(SigmaX) <- diag(SigmaY) <- 1

    
    for (k in 1:nsims) {
    
      cat("n =", n, "rho =", rho, "nsim =", k, "\n")
      seednum <-(i - 1) * J * K + (j - 1) * K + k
      writeLines(paste("i =", i, "j =", j, "k =", k),
                 con = "null_sim_log.txt")
      cat("------", seednum, "\n")
      set.seed(seednum)
      X <- mvrnormR(n, rep(0, dx), SigmaX)
      Y <- mvrnormR(n, rep(0, dy), SigmaY)
      
      set.seed(12345)
      chisq_res <- cbce(X, Y, backend="chisq", exhaustive=TRUE,
                        init_method="no-multiple-testing")
      
      set.seed(12345)
      indepChiSq_res <- cbce(X, Y, backend="indepChiSq", exhaustive=TRUE,
                             init_method="no-multiple-testing")
      
      # Extracting metrics
      stableIndxs1 <- indepChiSq_res$finalIndxs
      stableIndxs2 <- chisq_res$finalIndxs
      meth1_ncomm_k <- length(stableIndxs1)
      meth2_ncomm_k <- length(stableIndxs2)
      meth1_nnode_k <- length(unlist(lapply(stableIndxs1, function (i) {
        indepChiSq_res$extract_res[[i]]$StableComm
      })))
      meth2_nnode_k <- length(unlist(lapply(stableIndxs1, function (i) {
        chisq_res$extract_res[[i]]$StableComm
      })))
      
      # Storing metrics
      meth1_ncomm[i, j] <- meth1_ncomm[i, j] + meth1_ncomm_k / nsims
      meth2_ncomm[i, j] <- meth2_ncomm[i, j] + meth2_ncomm_k / nsims
      meth1_nnode[i, j] <- meth1_nnode[i, j] + meth1_nnode_k / nsims
      meth2_nnode[i, j] <- meth2_nnode[i, j] + meth2_nnode_k / nsims
      diffs_ncomm[i, j] <- diffs_ncomm[i, j] + (meth2_ncomm_k - meth1_ncomm_k) / nsims
      diffs_nnode[i, j] <- diffs_nnode[i, j] + (meth2_nnode_k - meth1_nnode_k) / nsims      
      
    }
    
  }
  
}

save(meth1_ncomm, meth2_ncomm, meth1_nnode, meth2_nnode,
     diffs_ncomm, diffs_nnode, file="null_sims_data.RData")

