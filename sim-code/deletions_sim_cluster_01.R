##### Runs model on simulated data (for MAF interval (0.1, 0.105))

args <- commandArgs(trailingOnly=TRUE)
iter.number <- as.integer(args[1])

set.seed(25)
library(dplyr)
library(purrr)
library(parallel)
source("deletion_sim_functions.R")
source("deletion_lik_functions.R")
n.cores <- detectCores()

### Resample trios to obtain bootstrapped trio genotype counts
get_boot_G <- function(G.df) {
  sampled_trios <- sample(1:nrow(G.df), replace = TRUE)
  G.list <- colSums(G.df[sampled_trios,])
  return(G.list)
}

### Helper function used in optim. Negative log-likelihood using the model with or without deletions is calculated.
### When model = "null", theta is a vector of six parameters. When model = "alt", theta is a vector of seven parameters.    
optim_function <- function(theta, Pi_parents, G, model) {
  if (model == "reduced.null") {
    res <- neg_sum_log_lik(exp(theta[1]), exp(theta[1])/100, exp(theta[2]), exp(theta[3]), 
                           exp(theta[4])/100, exp(theta[4]), 0, 0,
                           0, 0, 0, Pi_parents, G)
  }
  if (model == "null") {
    res <- neg_sum_log_lik(exp(theta[1]), exp(theta[2]), exp(theta[3]), exp(theta[4]), 
                           exp(theta[5]), exp(theta[6]), 0, 0,
                           0, 0, 0, Pi_parents, G)
  }
  if (model == "reduced.alt") {
    res <- neg_sum_log_lik(exp(theta[1]), exp(theta[1])/100, exp(theta[2]), exp(theta[3]), 
                           exp(theta[4])/100, exp(theta[4]), exp(theta[1]), exp(theta[1])/100,
                           exp(theta[4])/100, exp(theta[4]), exp(theta[5]), Pi_parents, G)
  }
  if (model == "alt") {
    res <- neg_sum_log_lik(exp(theta[1]), exp(theta[2]), exp(theta[3]), exp(theta[4]), 
                           exp(theta[5]), exp(theta[6]), exp(theta[1]), exp(theta[2]),
                           exp(theta[5]), exp(theta[6]), exp(theta[7]), Pi_parents, G)
  }
  if (model == "full.alt") {
    res <- neg_sum_log_lik(exp(theta[1]), exp(theta[2]), exp(theta[3]), exp(theta[4]), 
                           exp(theta[5]), exp(theta[6]), exp(theta[7]), exp(theta[8]),
                           exp(theta[9]), exp(theta[10]), exp(theta[11]), Pi_parents, G)
  }
  return(res)
}

### Returns the estimated ancestral parental genotype frequencies and optimization results for a single MAF interval
### G: vector of trio genotype counts for a single MAF interval
### Pi_parents_anc: Ancestral parental genotype frequencies (not used)
### par: specifies whether to return just the estimated parameter values, or the entire optim output
run_real_trios <- function(G, Pi_parents_anc, par = FALSE) {
  G <- as.numeric(G)
  Pi_parents <- get_Pi_parents_consistent(G)
  
  start <- Sys.time()
  print("starting optim")
  
  res_optim10000.null <- optim(c(rep(log(5e-4), 6)), optim_function,
                               Pi_parents = Pi_parents, G = G, model = "null",
                               method = "SANN", control = list(maxit = 10000))
  res_optim10000 <- optim(c(rep(log(5e-4), 7)), optim_function,
                          Pi_parents = Pi_parents, G = G, model = "alt",
                          method = "SANN", control = list(maxit = 10000))
  
  print("time spent on optim: ")
  print( Sys.time() - start )
  
  if (par == TRUE) {return(c(res_optim10000.null$par, res_optim10000$par))}
  return(list(res_optim10000.null, res_optim10000))
}

### Simulates trio genotypes and returns point estimates and bootstrapped estimates
### n_boot: number of bootstrap iterations
### p: start of MAF interval (p, p + 0.005)
do_one <- function(n_boot, n_trios, n_pos, n_reps, p, theta, gamma) {
  start_one <- Sys.time()
  ### Simulate trio genotype counts
  d_ <- gen_dataset(n_trios, n_pos, n_reps, p, theta, gamma)
  
  ### Obtain raw and bootstrapped trio genotype counts
  G.df <- d_[[1]]
  obs_G <- colSums(G.df)
  boot_G <- replicate(n_boot, get_boot_G(G.df))
  boot_G_lst <- lapply(seq_len(ncol(boot_G)), function(i) boot_G[,i])
  
  Pi_parents_anc <- d_[[2]]
  
  ### Fit model to raw and bootstrapped trio genotype counts
  res.point <- run_real_trios(obs_G, Pi_parents_anc)
  res.boot <- lapply(boot_G_lst, run_real_trios,
                     Pi_parents_anc = Pi_parents_anc, 
                     par = TRUE) %>% unlist() %>% matrix(byrow = TRUE, nrow = n_boot)
  
  print("time spent on iteration:")
  print( Sys.time() - start_one )
  return(list(res.point, res.boot, G.df, Pi_parents_anc))
}

theta01 <- 2e-4
theta02 <- 5e-7
theta00 <- 1 - theta01 - theta02

theta10 <- 9e-5
theta12 <- 2e-4
theta11 <- 1 - theta10 - theta12

theta20 <- 2e-6
theta21 <- 9e-4
theta22 <- 1 - theta20 - theta21

alpha <- 1

theta31 <- alpha*theta01
theta32 <- alpha*theta02
theta30 <- 1 - theta31 - theta32

theta40 <- alpha*theta20
theta41 <- alpha*theta21
theta42 <- 1 - theta40 - theta41

theta <- c(theta00, theta01, theta02, 
           theta10, theta11, theta12, 
           theta20, theta21, theta22,
           theta30, theta31, theta32,
           theta40, theta41, theta42)

gamma <- 2e-4

n_trios <- 100 
n_pos <- 10^5
n_reps <- 3000

num.iter.vec <- rep(100, n.cores)

options(scipen = 999)

file.name0.01 <- paste0("res.obsPi.boot.", paste(theta01, theta02, theta10, theta12, theta20, theta21, gamma, sep = '.'),
                        ".p.", 0.01,
                        ".n_trios.", n_trios,
                        ".n_pos.", n_pos*n_reps,
                        ".iter_optim", 10000, ".",
                        iter.number, ".rds")
res0.01 <- mclapply(num.iter.vec, do_one, 
                    n_trios = n_trios, 
                    n_pos = n_pos, 
                    n_reps = n_reps,
                    p = 0.01, 
                    theta = theta, 
                    gamma = gamma,
                    mc.cores = n.cores)
saveRDS(res0.01, file.name0.01)