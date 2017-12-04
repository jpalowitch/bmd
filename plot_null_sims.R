#Args <- commandArgs(TRUE)
Args <- c("n")
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

# ------------------------------------------------------------------------------

load(paste0("null_sims_data_", Type, ".RData"))

max_comms = max(c(meth1_ncomm, meth2_ncomm))
max_nodes = max(c(meth1_nnode, meth2_nnode))

rownames(meth1_ncomm) <- rownames(meth2_ncomm) <- pars
rownames(meth1_nnode) <- rownames(meth2_nnode) <- pars
colnames(meth1_ncomm) <- colnames(meth2_ncomm) <- rhos
colnames(meth1_nnode) <- colnames(meth2_nnode) <- rhos

library(ggplot2)
library(reshape2)

nL <- length(ns); nR <- length(rhos)
methvec <- c(rep("strong_chisq", nL * nR), rep("weak_chisq", nL * nR)) 
ncomm_df <- rbind(melt(meth1_ncomm), melt(meth2_ncomm))
ncomm_df$Method <- methvec
nnode_df <- rbind(melt(meth1_nnode), melt(meth2_nnode))
nnode_df$Method <- methvec
colnames(nnode_df)[1:3] <- c("par", "rho", "Avg.Num.Nodes")
colnames(ncomm_df)[1:3] <- c("par", "rho", "Avg.Num.Comms")

head(ncomm_df)

p1 <- ggplot(ncomm_df, aes(x=par, y=rho, fill=Avg.Num.Comms)) + 
  facet_wrap(~Method, ncol=2) + geom_tile() + xlab(Type) + ylab("rho") + 
  ggtitle("average number of communities")

p2 <- ggplot(nnode_df, aes(x=par, y=rho, fill=Avg.Num.Nodes)) + 
  facet_wrap(~Method, ncol=2) + geom_tile() + xlab(Type) + ylab("rho") + 
  ggtitle("average number of nodes")

ggsave(paste0("null_sims_plot_", Type, "_comms.png"), p1)
ggsave(paste0("null_sims_plot_", Type, "_nodes.png"), p2)
