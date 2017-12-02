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

load("null_sims_data.RData")

max_comms = max(c(meth1_ncomm, meth2_ncomm))
max_nodes = max(c(meth1_nnode, meth2_nnode))

rownames(meth1_ncomm) <- rownames(meth2_ncomm) <- ns
rownames(meth1_nnode) <- rownames(meth2_nnode) <- ns
colnames(meth1_ncomm) <- rownames(meth2_ncomm) <- rhos
colnames(meth1_nnode) <- rownames(meth2_nnode) <- rhos

library(ggplot2)
library(reshape2)

nL <- length(ns); nR <- length(rhos)
methvec <- c(rep("meth1", nL * nR), rep("meth2", nL * nR)) 
ncomm_df <- rbind(melt(meth1_ncomm), melt(meth2_ncomm))
ncomm_df$Method <- methvec
nnode_df <- rbind(melt(meth1_nnode), melt(meth2_nnode))
nnode_df$Method <- methvec
colnames(nnode_df)[1:3] <- c("n", "rho", "Avg.Num.Comm.Nodes")
colnames(ncomm_df)[1:3] <- c("n", "rho", "Avg.Num.Comms")

head(ncomm_df)

p <- ggplot(ncomm_df, aes(x = n, y = rho, fill = Avg.Num.Comms)) + 
  facet_wrap(~Method, ncol = 2) + geom_tile()
