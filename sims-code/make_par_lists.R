source("sim_eQTL_network.R")

# For the experiments, do you want to shove the dec vec up to one?
shove_dec <- TRUE

# Global settings
par_divs <- 9
par_seq_dec <- 1:par_divs / (par_divs + 1)
par_seq  <- round(100 * (1:par_divs / (par_divs + 1)))
par_dirs <- as.character(par_seq)

# Give the names of your experiments: must be manually entered.
total_expers <- as.character(1:17)

if (!dir.exists("sims-results/sbm-par-lists"))
  dir.create("sims-results/sbm-par-lists", recursive = TRUE)

# Experiment 1 -----------------------------------------------------------------

main_text <- "Increase n"
par_list <- make_param_list()
pars <- c("n")
xlab <- "Sample Size"
palpha <- FALSE
axis_par <- 1
par_settings <- matrix(0, 1, par_divs)
par_settings[1, ] <- round(200 * (par_seq_dec + min(par_seq_dec) * 
                                     as.numeric(shove_dec)))

save(par_list,
     main_text,
     axis_par,
     pars,
     xlab, palpha,
     par_settings,
     par_seq,
     par_divs,
     par_dirs,
     file = "sims-results/sbm-par-lists/experiment1.RData")

# Experiment 2 -----------------------------------------------------------------

main_text <- "Decrease min bimodule-half size"
par_list <- make_param_list()
pars <- c("cmin")
xlab <- "Min BM-half size"
palpha <- FALSE
axis_par <- 1
par_settings <- matrix(0, 1, par_divs)
par_settings[1, ] <- round(50 * (par_seq_dec + min(par_seq_dec) * 
                                    as.numeric(shove_dec)))

save(par_list,
     main_text,
     axis_par,
     pars,
     xlab, palpha,
     par_settings,
     par_seq,
     par_divs,
     par_dirs,
     file = "sims-results/sbm-par-lists/experiment2.RData")

# Experiment 3 -----------------------------------------------------------------

main_text <- "Increase amount of bg"
par_list <- make_param_list()
pars <- c("bgmult")
xlab <- "#BG vars/#BM vars"
palpha <- FALSE
axis_par <- 1
par_settings <- matrix(0, 1, par_divs)
par_settings[1, ] <- c(0, 0.5, 10 * (par_seq_dec[-c(1, par_divs)]) - 1)

save(par_list,
     main_text,
     axis_par,
     pars,
     xlab, palpha,
     par_settings,
     par_seq,
     par_divs,
     par_dirs,
     file = "sims-results/sbm-par-lists/experiment3.RData")

# Experiment 4 -----------------------------------------------------------------

main_text <- "Decrease mean of beta parameters"
par_list <- make_param_list()
pars <- c("cormean")
xlab <- "Mean of beta params"
palpha <- FALSE
axis_par <- 1
par_settings <- matrix(0, 1, par_divs)
par_settings[1, ] <- 1 * (par_seq_dec + min(par_seq_dec) * as.numeric(shove_dec))

save(par_list,
     main_text,
     axis_par,
     pars,
     xlab, palpha,
     par_settings,
     par_seq,
     par_divs,
     par_dirs,
     file = "sims-results/sbm-par-lists/experiment4.RData")

# Experiment 5 -----------------------------------------------------------------

main_text <- "Decrease eQTL probability"
par_list <- make_param_list()
pars <- c("p")
xlab <- "eQTL probability"
palpha <- FALSE
axis_par <- 1
par_settings <- matrix(0, 1, par_divs)
par_settings[1, ] <- 0.5 * (par_seq_dec + min(par_seq_dec) * as.numeric(shove_dec))

save(par_list,
     main_text,
     axis_par,
     pars,
     xlab, palpha,
     par_settings,
     par_seq,
     par_divs,
     par_dirs,
     file = "sims-results/sbm-par-lists/experiment5.RData")

# Experiment 6 -----------------------------------------------------------------

main_text <- "Increase intra-X correlation"
par_list <- make_param_list()
pars <- c("rho")
xlab <- "Intra-X corr"
palpha <- FALSE
axis_par <- 1
par_settings <- matrix(0, 1, par_divs)
par_settings[1, ] <- 0.3 * (par_seq_dec) - 0.03

save(par_list,
     main_text,
     axis_par,
     pars,
     xlab, palpha,
     par_settings,
     par_seq,
     par_divs,
     par_dirs,
     file = "sims-results/sbm-par-lists/experiment6.RData")

# Experiment 7 -----------------------------------------------------------------

main_text <- "Increase base noise scaling"
par_list <- make_param_list()
pars <- c("s2")
xlab <- "Noise scaling"
palpha <- FALSE
axis_par <- 1
par_settings <- matrix(0, 1, par_divs)
par_settings[1, ] <- 50 * (par_seq_dec + min(par_seq_dec) * as.numeric(shove_dec)) - 5

save(par_list,
     main_text,
     axis_par,
     pars,
     xlab, palpha,
     par_settings,
     par_seq,
     par_divs,
     par_dirs,
     file = "sims-results/sbm-par-lists/experiment7.RData")

# Experiment 8 -----------------------------------------------------------------

main_text <- "Increase amount of bg, w/corNoiseX"
par_list <- make_param_list(corNoiseX = TRUE)
pars <- c("bgmult")
xlab <- "#BG vars/#BM vars"
palpha <- FALSE
axis_par <- 1
par_settings <- matrix(0, 1, par_divs)
par_settings[1, ] <- c(0, 0.5, 10 * (par_seq_dec[-c(1, par_divs)]) - 1)

save(par_list,
     main_text,
     axis_par,
     pars,
     xlab, palpha,
     par_settings,
     par_seq,
     par_divs,
     par_dirs,
     file = "sims-results/sbm-par-lists/experiment8.RData")

# Experiment 9 -----------------------------------------------------------------

main_text <- "Increase amount of bg, w/corNoiseY"
par_list <- make_param_list(corNoiseY = TRUE)
pars <- c("bgmult")
xlab <- "#BG vars/#BM vars"
palpha <- FALSE
axis_par <- 1
par_settings <- matrix(0, 1, par_divs)
par_settings[1, ] <- c(0, 0.5, 10 * (par_seq_dec[-c(1, par_divs)]) - 1)

save(par_list,
     main_text,
     axis_par,
     pars,
     xlab, palpha,
     par_settings,
     par_seq,
     par_divs,
     par_dirs,
     file = "sims-results/sbm-par-lists/experiment9.RData")

# Experiment 10 -----------------------------------------------------------------

main_text <- "Increase #bg w/corNoiseX&Y"
par_list <- make_param_list(corNoiseX = TRUE, corNoiseY = TRUE)
pars <- c("bgmult")
xlab <- "#BG vars/#BM vars"
palpha <- FALSE
axis_par <- 1
par_settings <- matrix(0, 1, par_divs)
par_settings[1, ] <- c(0, 0.5, 10 * (par_seq_dec[-c(1, par_divs)]) - 1)

save(par_list,
     main_text,
     axis_par,
     pars,
     xlab, palpha,
     par_settings,
     par_seq,
     par_divs,
     par_dirs,
     file = "sims-results/sbm-par-lists/experiment10.RData")

# Experiment 11 -----------------------------------------------------------------
#----** Special experiment which adjusts algorithm parameter **----#

main_text <- "Decrease alpha"
par_list <- make_param_list()
pars <- c("alpha")
xlab <- expression(alpha)
palpha <- TRUE
axis_par <- 1
par_settings <- matrix(0, 1, par_divs)
par_settings[1, ] <- 0.2 * 2^(-c(0:(par_divs - 1)))

save(par_list,
     main_text,
     axis_par,
     pars,
     xlab, palpha,
     par_settings,
     par_seq,
     par_divs,
     par_dirs,
     file = "sims-results/sbm-par-lists/experiment11.RData")

# Experiment 12 -----------------------------------------------------------------

main_text <- "Increase n"
par_list <- make_param_list()
pars <- c("n")
xlab <- "Sample Size"
palpha <- FALSE
axis_par <- 1
par_settings <- matrix(0, 1, par_divs)
par_settings[1, ] <- round(2000 * (par_seq_dec + min(par_seq_dec) * 
                                     as.numeric(shove_dec)))

save(par_list,
     main_text,
     axis_par,
     pars,
     xlab, palpha,
     par_settings,
     par_seq,
     par_divs,
     par_dirs,
     file = "sims-results/sbm-par-lists/experiment12.RData")

# Experiment 13 -----------------------------------------------------------------

main_text <- "Increase n (w/ clearer signal)"
par_list <- make_param_list(p = 1, rho = 0.5, cormean = 3, s2 = 1, bgmult = 0,
                            cmin = 100, b = 4)
pars <- c("n")
xlab <- "Sample Size"
palpha <- FALSE
axis_par <- 1
par_settings <- matrix(0, 1, par_divs)
par_settings[1, ] <- round(3000 * (par_seq_dec + min(par_seq_dec) * 
                                     as.numeric(shove_dec)))

save(par_list,
     main_text,
     axis_par,
     pars,
     xlab, palpha,
     par_settings,
     par_seq,
     par_divs,
     par_dirs,
     file = "sims-results/sbm-par-lists/experiment13.RData")

# Experiment 14 -----------------------------------------------------------------
#----** Special experiment which adjusts algorithm parameter **----#

main_text <- "Decrease alpha (w/ clearer signal)"
par_list <- make_param_list(p = 1, rho = 0.5, cormean = 3, s2 = 1, bgmult = 0,
                            cmin = 100, b = 4, n = 1000)
pars <- c("alpha")
xlab <- expression(alpha)
palpha <- TRUE
axis_par <- 1
par_settings <- matrix(0, 1, par_divs)
par_settings[1, ] <- 0.2 * 2^(-c(0:(par_divs - 1)))

save(par_list,
     main_text,
     axis_par,
     pars,
     xlab, palpha,
     par_settings,
     par_seq,
     par_divs,
     par_dirs,
     file = "sims-results/sbm-par-lists/experiment14.RData")

# Experiment 15 -----------------------------------------------------------------

main_text <- "Increase n, bm-size power law"
par_list <- make_param_list(size_law = -1, cmin = 10, cmax = 100)
pars <- c("n")
xlab <- "Sample Size"
palpha <- FALSE
axis_par <- 1
par_settings <- matrix(0, 1, par_divs)
par_settings[1, ] <- round(200 * (par_seq_dec + min(par_seq_dec) * 
                                    as.numeric(shove_dec)))

save(par_list,
     main_text,
     axis_par,
     pars,
     xlab, palpha,
     par_settings,
     par_seq,
     par_divs,
     par_dirs,
     file = "sims-results/sbm-par-lists/experiment15.RData")

# Experiment 16 -----------------------------------------------------------------

main_text <- "Increase variance, bm-size power law"
par_list <- make_param_list(size_law = -1, cmin = 10, cmax = 100)
pars <- c("s2")
xlab <- "Decrease mean of beta params, bm-size power law"
palpha <- FALSE
axis_par <- 1
par_settings <- matrix(0, 1, par_divs)
par_settings[1, ] <- 20 * (par_seq_dec + min(par_seq_dec) * as.numeric(shove_dec))

save(par_list,
     main_text,
     axis_par,
     pars,
     xlab, palpha,
     par_settings,
     par_seq,
     par_divs,
     par_dirs,
     file = "sims-results/sbm-par-lists/experiment16.RData")

# Experiment 17 -----------------------------------------------------------------

main_text <- "Decrease eQTL probability, bm-size power law"
par_list <- make_param_list(size_law = -1, cmin = 10, cmax = 100)
pars <- c("p")
xlab <- "eQTL probability"
palpha <- FALSE
axis_par <- 1
par_settings <- matrix(0, 1, par_divs)
par_settings[1, ] <- 0.5 * (par_seq_dec + min(par_seq_dec) * as.numeric(shove_dec))

save(par_list,
     main_text,
     axis_par,
     pars,
     xlab, palpha,
     par_settings,
     par_seq,
     par_divs,
     par_dirs,
     file = "sims-results/sbm-par-lists/experiment17.RData")


writeLines(total_expers, "sims-results/exper-names.txt")


writeLines(total_expers, "sims-results/exper-names.txt")


