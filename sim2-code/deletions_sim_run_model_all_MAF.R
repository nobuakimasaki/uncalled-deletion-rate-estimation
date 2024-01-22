##### Runs model on simulated trio genotype counts from 100 MAF intervals 

code_start_time <- Sys.time()

args <- commandArgs(trailingOnly=TRUE)
boot.iter <- args[1]
iter.number <- args[2]

print("num bootstrap iterations: ")
print(boot.iter)

print("iter number: ")
print(iter.number)

set.seed(25)
library(dplyr)
library(purrr)
library(parallel)
source("deletion_sim_functions.R")
source("deletion_lik_functions.R")
n.cores <- detectCores()

### Helper function used in optim. Negative log-likelihood using the model with or without deletions is calculated.
### When model = "null", theta is a vector of six parameters. When model = "alt", theta is a vector of seven parameters.    
optim_function <- function(theta, Pi_parents, G, model) {
  if (model == "null") {
    res <- neg_sum_log_lik(exp(theta[1]), exp(theta[2]), exp(theta[3]), exp(theta[4]), 
                           exp(theta[5]), exp(theta[6]), 0, 0,
                           0, 0, 0, Pi_parents, G)
  }
  if (model == "alt") {
    res <- neg_sum_log_lik(exp(theta[1]), exp(theta[2]), exp(theta[3]), exp(theta[4]), 
                           exp(theta[5]), exp(theta[6]), exp(theta[1]), exp(theta[2]),
                           exp(theta[5]), exp(theta[6]), exp(theta[7]), Pi_parents, G)
  }
  return(res)
}

### Returns the estimated ancestral parental genotype frequencies and optimization results for a single MAF interval
### G: vector of trio genotype counts for a single MAF interval
run_real_trios <- function(G) {
  G <- as.numeric(G)
  Pi_parents <- get_Pi_parents_consistent(G)
  
  start <- Sys.time()
  print("starting optim")
  
  # res_optim10000.1 <- optim(c(rep(log(5e-4), 7)), optim_function,
  #                         Pi_parents = Pi_parents, G = G, model = "alt",
  #                         method = "Nelder-Mead")
  # res_optim10000.2 <- optim(c(rep(log(5e-4), 7)), optim_function,
  #                         Pi_parents = Pi_parents, G = G, model = "alt",
  #                         method = "BFGS")
  # res_optim10000.3 <- optim(c(rep(log(5e-4), 7)), optim_function,
  #                         Pi_parents = Pi_parents, G = G, model = "alt",
  #                         method = "CG")
  # res_optim10000.4 <- optim(c(rep(log(5e-4), 7)), optim_function,
  #                         Pi_parents = Pi_parents, G = G, model = "alt",
  #                         method = "L-BFGS-B")
  res_optim10000 <- optim(c(rep(log(5e-4), 7)), optim_function,
                          Pi_parents = Pi_parents, G = G, model = "alt",
                          method = "SANN", control = list(maxit = 10000))
  
  # res_optim10000.null.1 <- optim(c(rep(log(5e-4), 6)), optim_function,
  #                           Pi_parents = Pi_parents, G = G, model = "null",
  #                           method = "Nelder-Mead")
  # res_optim10000.null.2 <- optim(c(rep(log(5e-4), 6)), optim_function,
  #                           Pi_parents = Pi_parents, G = G, model = "null",
  #                           method = "BFGS")
  # res_optim10000.null.3 <- optim(c(rep(log(5e-4), 6)), optim_function,
  #                           Pi_parents = Pi_parents, G = G, model = "null",
  #                           method = "CG")
  # res_optim10000.null.4 <- optim(c(rep(log(5e-4), 6)), optim_function,
  #                           Pi_parents = Pi_parents, G = G, model = "null",
  #                           method = "L-BFGS-B")
  # res_optim10000.null <- optim(c(rep(log(5e-4), 6)), optim_function,
  #                              Pi_parents = Pi_parents, G = G, model = "null",
  #                              method = "SANN", control = list(maxit = 10000))
  
  print("time spent on optim: ")
  print( Sys.time() - start )
  
  return(list(Pi_parents,
              res_optim10000))
}

### Fits the model to all MAF intervals. 
### trios: In the case of a bootstrap sample, trios is the vector of resampled trio ids. 
### If we're fitting the model to the original sample, trios is the vector of unique trio ids in the dataset.
boot_one <- function(trios, G.df) {
  ### We accumulate all of the trio genotypes from each trio id
  trio_counts <- matrix(nrow = 0, ncol = 19) %>% as.data.frame()
  colnames(trio_counts) <- c(paste0("V", seq(1,18)), "MAF_start")
  for (tr in trios) {
    temp <- G.df %>% filter(trio_id == tr) %>% select(-trio_id)
    trio_counts <- rbind(trio_counts, temp)
  }
  trio_counts <- trio_counts %>% group_by(MAF_start) %>% summarise_all(sum) %>% select(-MAF_start, V1, V2, V3, V4, V5, V6,
                                                                                       V7, V8, V9, V10, V11, V12,
                                                                                       V13, V14, V15, V16, V17, V18) 
  
  ### G.list is split according to the MAF intervals
  G.list <- split(trio_counts, seq(nrow(trio_counts)))
  
  start.lapply <- Sys.time()
  res <- mclapply(G.list, run_real_trios, mc.cores = 24)
  print("total time spent on optim (on all maf intervals): ")
  print( Sys.time() - start.lapply )
  
  return(res)
}

### Simulated observed trio genotype counts
G.df <- readRDS("res\\sim_obs_trio_cnts_all_MAF.rds")
# c_parents_anc <- readRDS("sim_predel_parental_cnts_all_MAF.rds")
# c_trios_true <- readRDS("sim_true_trio_cnts_all_MAF.rds")

### Resample trios and run the model on all MAF intervals in the resulting dataset
unique.trios <- unique(G.df$trio_id)
trios.resamp <- replicate(as.numeric(boot.iter) + 1, sample(unique.trios, length(unique.trios), replace = TRUE)) %>%
  t()
trios.resamp.list <- split(trios.resamp, seq(nrow(trios.resamp)))
trios.resamp.list[[1]] <- unique.trios

res <- mclapply(trios.resamp.list, boot_one, 
                G.df = G.df, 
                mc.cores = 1)

resampled.trios.file.name <- paste0("resamp.trios.",
                                    boot.iter, ".",
                                    iter.number, ".all_MAF.sim.rds")
est.res.file.name <- paste0("est.res.",
                            boot.iter, ".",
                            iter.number, ".all_MAF.sim.rds")

saveRDS(trios.resamp.list, resampled.trios.file.name)
saveRDS(res, est.res.file.name)

code_end_time <- Sys.time()
print("time elapsed: ")
print(code_end_time - code_start_time)