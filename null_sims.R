 Args <- commandArgs(TRUE)
#Args <- c("n")
if (length(Args) < 1) {
  stop('usage: [type]\n')
}

Type = Args[1]

if (!Type %in% c("n", "bg")) {
  stop('usage: [type]\n[type] must be in c("n", "bg")')
}

# ------------------------------------------------------------------------------

source("mvrnormR.R")
library(cbce)
source("null_sims_control.R")

# ------------------------------------------------------------------------------

# Setting up 
if (Type == "n") {
  pars = ns
} else if (Type == "bg") {
  pars = bgs
}
meth1_ncomm <- meth2_ncomm <- diffs_ncomm <- matrix(0, length(pars), length(rhos))
meth1_nnode <- meth2_nnode <- diffs_nnode <- matrix(0, length(pars), length(rhos))
I <- length(pars)
J <- length(rhos)
K <- nsims

# ------------------------------------------------------------------------------

for (i in seq_along(pars)) {
  
  n <- ns[i]
  
  if (Type == "n") {
    bg <- bgs[1]
    n <- pars[i]
  } else if (Type == "bg") {
    bg <- pars[i]
    n <- ns[2]
  }
  dx <- mx + bg; dy <- my + bg
  SigmaX <- matrix(0, dx, dx)
  SigmaY <- matrix(0, dy, dy)
  
  for (j in seq_along(rhos)) {
    
    rho <- rhos[j]
    SigmaX[1:mx, 1:mx] <- SigmaY[1:mx, 1:mx] <- rho
    diag(SigmaX) <- diag(SigmaY) <- 1

    
    for (k in 1:nsims) {
    
      if (Type == "n") {
        cat("n =", n, "rho =", rho, "nsim =", k, "\n")
      } else if (Type == "bg") {
        cat("bg =", bg, "rho =", rho, "nsim =", k, "\n")
      }
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
     diffs_ncomm, diffs_nnode, file=paste0("null_sims_data_", Type, ".RData"))

