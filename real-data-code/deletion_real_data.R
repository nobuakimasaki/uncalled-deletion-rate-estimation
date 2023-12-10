##### Code to run model on the UK Biobank sequence data.

code_start_time <- Sys.time()

### Argument is simply the UK Biobank data file name on the cluster
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("Arguments must be supplied", call.=FALSE) }

input.file <- args[1]

print("input file: ")
print(input.file)

set.seed(25)
source("deletion_lik_functions.R")
library(dplyr)
library(purrr)
library(parallel)
n.cores <- detectCores()

### Helper functions to add matrices and undo cumulative sum
add <- function(x) Reduce("+", x)

undo_cumsum <- function(x) {
  c(x[1],diff(x))
}

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
  
  res_optim10000.1 <- optim(c(rep(log(5e-4), 7)), optim_function,
                            Pi_parents = Pi_parents, G = G, model = "alt",
                            method = "Nelder-Mead")
  res_optim10000.2 <- optim(c(rep(log(5e-4), 7)), optim_function,
                            Pi_parents = Pi_parents, G = G, model = "alt",
                            method = "BFGS")
  res_optim10000.3 <- optim(c(rep(log(5e-4), 7)), optim_function,
                            Pi_parents = Pi_parents, G = G, model = "alt",
                            method = "CG")
  res_optim10000.4 <- optim(c(rep(log(5e-4), 7)), optim_function,
                            Pi_parents = Pi_parents, G = G, model = "alt",
                            method = "L-BFGS-B")
  res_optim10000.5 <- optim(c(rep(log(5e-4), 7)), optim_function,
                            Pi_parents = Pi_parents, G = G, model = "alt",
                            method = "SANN", control = list(maxit = 10000))
  
  res_optim10000.null.1 <- optim(c(rep(log(5e-4), 6)), optim_function,
                                 Pi_parents = Pi_parents, G = G, model = "null",
                                 method = "Nelder-Mead")
  res_optim10000.null.2 <- optim(c(rep(log(5e-4), 6)), optim_function,
                                 Pi_parents = Pi_parents, G = G, model = "null",
                                 method = "BFGS")
  res_optim10000.null.3 <- optim(c(rep(log(5e-4), 6)), optim_function,
                                 Pi_parents = Pi_parents, G = G, model = "null",
                                 method = "CG")
  res_optim10000.null.4 <- optim(c(rep(log(5e-4), 6)), optim_function,
                                 Pi_parents = Pi_parents, G = G, model = "null",
                                 method = "L-BFGS-B")
  res_optim10000.null.5 <- optim(c(rep(log(5e-4), 6)), optim_function,
                                 Pi_parents = Pi_parents, G = G, model = "null",
                                 method = "SANN", control = list(maxit = 10000))
  
  print("time spent on optim: ")
  print( Sys.time() - start )
  
  return(list(Pi_parents,
              res_optim10000.1, res_optim10000.2, res_optim10000.3, res_optim10000.4, res_optim10000.5,
              res_optim10000.null.1, res_optim10000.null.2, res_optim10000.null.3, res_optim10000.null.4, res_optim10000.null.5))
}

### Fits the model to all MAF intervals. 
### trios: In the case of a bootstrap sample, trios is the vector of resampled trio ids. 
### If we're fitting the model to the original sample, trios is the vector of unique trio ids in the dataset.
boot_one <- function(trios, G.df) {
  n_intervals <- G.df %>% filter(V2 == trios[1]) %>% nrow() 
  n_intervals <- n_intervals - 1
  
  ### We accumulate all of the trio genotypes from each trio id
  total <- rep(0, n_intervals*18) %>% matrix(nrow = n_intervals) %>% as.data.frame()
  for (tr in trios) {
    G.filt <- G.df %>% filter(V2 == tr) %>% arrange(V1)
    values <- G.filt[,3:20] %>% 
      mutate_all(function(x) as.numeric(as.character(x))) %>%
      apply(2, undo_cumsum) %>% 
      as.data.frame() %>%
      mutate_all(function(x) as.numeric(as.character(x)))
    ### We remove the first MAF interval which is (0, 0]
    values <- values[-1,]
    total <- add(list(total, values))
  }
  
  ### G.list is split according to the MAF intervals
  G.list <- split(total, seq(nrow(total)))
  
  start.lapply <- Sys.time()
  res <- mclapply(G.list, run_real_trios, mc.cores = 24)
  print("total time spent on optim (on all maf intervals): ")
  print( Sys.time() - start.lapply )
  
  return(res)
}

### Load in the UK Biobank sequence data
ukbio <- read.table(input.file, header = TRUE)
G.df <- ukbio[,c(1,3,4:21)] %>% 
  unlist() %>% 
  matrix(nrow = 20, byrow = TRUE) %>% 
  t() %>%
  as.data.frame() %>%
  filter(V2 != "SUM" & V2 != "EXP" & V2 != "OBS") %>%
  arrange(V1)

est.res.file.name <- paste0("est.res.",
                            "pt", ".rds")

### Get unique trio ids
unique.trios <- unique(G.df$V2)
res <- boot_one(unique.trios, G.df)

saveRDS(res, est.res.file.name)

code_end_time <- Sys.time()
print("time elapsed: ")
print(code_end_time - code_start_time)
