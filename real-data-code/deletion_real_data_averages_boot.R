##### Code to obtain bootstrapped estimates for the average genotype error rate and deletion rate. 
##### Requires the file of bootstrapped estimates and resampled trios for each bootstrap to be in the same directory.

### Load in packages
library(dplyr)
library(tibble)
library(tidyr)
library(purrr)

### Helper functions to add matrices together and undo cumulative sum
add <- function(x) Reduce("+", x)

undo_cumsum <- function(x) {
  c(x[1],diff(x))
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("Input file must be supplied", call.=FALSE) }

input.file <- args[1]

### Function returning estimates of true trio frequencies
### Pi_parents: estimates of ancestral parental frequencies
get_Pi_c <- function(Pi_parents, gamma) {
  r <- 1-2*gamma
  c <- gamma
  Pi00 <- Pi_parents[[1]]
  Pi01 <- Pi_parents[[2]]
  Pi02 <- Pi_parents[[3]]
  Pi11 <- Pi_parents[[4]]
  Pi12 <- Pi_parents[[5]]
  Pi22 <- Pi_parents[[6]]
  
  Pic00 <- Pi00*r^2
  Pic01 <- Pi01*r^2
  Pic02 <- Pi02*r^2
  Pic03 <- Pi00*4*r*c + Pi01*r*c
  Pic04 <- Pi02*2*r*c + Pi01*r*c
  
  Pic11 <- Pi11*r^2
  Pic12 <- Pi12*r^2
  Pic13 <- Pi01*2*r*c + Pi11*2*r*c
  Pic14 <- Pi12*2*r*c + Pi11*2*r*c
  
  Pic22 <- Pi22*r^2
  Pic23 <- Pi02*2*r*c + Pi12*r*c
  Pic24 <- Pi22*4*r*c + Pi12*r*c
  
  Pic33 <- Pi00*4*c^2 + Pi01*2*c^2 + Pi11*c^2
  Pic34 <- Pi02*4*c^2 + Pi01*2*c^2 + Pi12*2*c^2 + Pi11*2*c^2
  
  Pic44 <- Pi22*4*c^2 + Pi12*2*c^2 + Pi11*c^2
  
  Pi000 = Pic00
  Pi010 = Pic01 / 2
  Pi011 = Pic01 / 2
  Pi021 = Pic02
  Pi030 = Pic03 / 2
  Pi033 = Pic03 / 2
  Pi041 = Pic04 / 2
  Pi043 = Pic04 / 2
  Pi110 = Pic11 / 4
  Pi111 = Pic11 / 2
  Pi112 = Pic11 / 4
  Pi121 = Pic12 / 2
  Pi122 = Pic12 / 2
  Pi130 = Pic13 / 4
  Pi131 = Pic13 / 4
  Pi133 = Pic13 / 4
  Pi134 = Pic13 / 4
  Pi141 = Pic14 / 4
  Pi142 = Pic14 / 4
  Pi143 = Pic14 / 4
  Pi144 = Pic14 / 4
  Pi222 = Pic22
  Pi231 = Pic23 / 2
  Pi234 = Pic23 / 2
  Pi242 = Pic24 / 2
  Pi244 = Pic24 / 2
  Pi330 = Pic33 / 4
  Pi333 = Pic33 / 2
  Pi335 = Pic33 / 4
  Pi341 = Pic34 / 4
  Pi343 = Pic34 / 4
  Pi344 = Pic34 / 4
  Pi345 = Pic34 / 4
  Pi442 = Pic44 / 4
  Pi444 = Pic44 / 2
  Pi445 = Pic44 / 4
  
  Pi0 <- matrix(c(Pi000, 0, 0, 0, 0,
                  Pi010, Pi011, 0, 0, 0,
                  0, Pi021, 0, 0, 0,
                  Pi030, 0, 0, Pi033, 0,
                  0, Pi041, 0, Pi043, 0), nrow = 5, ncol = 5, byrow = TRUE)/(1-Pi335-Pi345-Pi445)
  
  Pi1 <- matrix(c(0, 0, 0, 0, 0,
                  Pi110, Pi111, Pi112, 0, 0,
                  0, Pi121, Pi122, 0, 0,
                  Pi130, Pi131, 0, Pi133, Pi134,
                  0, Pi141, Pi142, Pi143, Pi144), nrow = 5, ncol = 5, byrow = TRUE)/(1-Pi335-Pi345-Pi445)
  
  Pi2 <- matrix(c(0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0,
                  0, 0, Pi222, 0, 0,
                  0, Pi231, 0, 0, Pi234,
                  0, 0, Pi242, 0, Pi244), nrow = 5, ncol = 5, byrow = TRUE)/(1-Pi335-Pi345-Pi445)
  
  Pi3 <- matrix(c(0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0,
                  Pi330, 0, 0, Pi333, 0,
                  0, Pi341, 0, Pi343, Pi344), nrow = 5, ncol = 5, byrow = TRUE)/(1-Pi335-Pi345-Pi445)
  
  Pi4 <- matrix(c(0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0,
                  0, 0, Pi442, 0, Pi444), nrow = 5, ncol = 5, byrow = TRUE)/(1-Pi335-Pi345-Pi445)
  
  Pi <- array(c(Pi0, Pi1, Pi2, Pi3, Pi4), dim = c(5,5,5))
  
  Pi <- aperm(Pi, c(3,1,2))
  return(Pi)
}

### Function returning error rate for a single MAF interval
### pnt: vector of point estimates
get_error_rate <- function(pnt, Pic) {
  Pik_vec <- rep(NA, 5)
  for (k in 1:5) {
    Pik <- 0
    for (i in 1:5) {
      for (j in i:5) {
        Pik <- Pik + Pic[i,j,k]
      }
    }
    Pik_vec[k] <- Pik
  }

  theta01 <- pnt[[1]]
  theta02 <- pnt[[2]]
  theta10 <- pnt[[3]]
  theta12 <- pnt[[4]]
  theta20 <- pnt[[5]]
  theta21 <- pnt[[6]]
  
  theta00 <- 1-theta01-theta02
  theta11 <- 1-theta10-theta12
  theta22 <- 1-theta20-theta21
  
  theta_vec <- c(theta00, theta11, theta22)
  
  delta <- 0
  for (k in 1:3) {
    delta <- delta + Pik_vec[k]*(1-theta_vec[k])
  }
  for (k in 4:5) {
    delta <- delta + Pik_vec[k]
  }
  return(delta)
}

### Load in bootstrapped estimates
res.1 <- readRDS("/projects/browning/brwnlab/masakin/deletions/Fall_2023/real_data_results/est.res.13.1.rds")
res.2 <- readRDS("/projects/browning/brwnlab/masakin/deletions/Fall_2023/real_data_results/est.res.13.2.rds")
res.3 <- readRDS("/projects/browning/brwnlab/masakin/deletions/Fall_2023/real_data_results/est.res.13.3.rds")
res.4 <- readRDS("/projects/browning/brwnlab/masakin/deletions/Fall_2023/real_data_results/est.res.13.4.rds")
res.5 <- readRDS("/projects/browning/brwnlab/masakin/deletions/Fall_2023/real_data_results/est.res.13.5.rds")
res.6 <- readRDS("/projects/browning/brwnlab/masakin/deletions/Fall_2023/real_data_results/est.res.13.6.rds")
res.7 <- readRDS("/projects/browning/brwnlab/masakin/deletions/Fall_2023/real_data_results/est.res.13.7.rds")
res.8 <- readRDS("/projects/browning/brwnlab/masakin/deletions/Fall_2023/real_data_results/est.res.13.8.rds")

# res.1 <- readRDS("result/est.res.13.1.rds")
# res.2 <- readRDS("result/est.res.13.2.rds")
# res.3 <- readRDS("result/est.res.13.3.rds")
# res.4 <- readRDS("result/est.res.13.4.rds")
# res.5 <- readRDS("result/est.res.13.5.rds")
# res.6 <- readRDS("result/est.res.13.6.rds")
# res.7 <- readRDS("result/est.res.13.7.rds")
# res.8 <- readRDS("result/est.res.13.8.rds")

### Each item in this large list represents a separate bootstrap iteration. The last four are removed because we're only using 100.
res.boot <- c(res.1,
              res.2,
              res.3,
              res.4,
              res.5,
              res.6,
              res.7,
              res.8)[-c(101,102,103,104)]

### Returns a matrix of point estimates, each row corresponding to a separate MAF interval
get_est <- function(pnt, i) {
  est.df <- lapply(pnt, function(x) {x[[i]][[1]]}) %>% unlist() %>% matrix(nrow = length(pnt), byrow = TRUE)  %>% exp() %>% as.data.frame()
  est.df <- est.df[-1,]
  return(est.df)
}

### Obtain matrix of point estimates for each bootstrap iteration
pnt.boot <- lapply(res.boot, get_est, i = 2)

### Returns a vector of gamma estimates, each entry corresponding to a separate MAF interval
get_gamma_list <- function(pnt, i) {
  est.df <- lapply(pnt, function(x) {x[[i]][[1]]}) %>% unlist() %>% matrix(nrow = length(pnt), byrow = TRUE)  %>% exp() %>% as.data.frame()
  gamma.list <- est.df[,7]
  gamma.list <- gamma.list[-1]
  return(gamma.list)
}

### Obtain vector of gamma estimates for each bootstrap iteration
gamma.boot <- lapply(res.boot, get_gamma_list, i = 2)

### Function returning a matrix of estimated ancestral parental frequencies, each row representing a separate MAF interval
get_Pi <- function(pnt) {
  est.df <- lapply(pnt, function(x) {x[[1]]}) %>% unlist() %>% matrix(nrow = length(pnt), byrow = TRUE) %>% as.data.frame()
  est.df <- est.df[-1,]
  return(est.df)
}

### Obtain matrix of estimated ancestral parental frequencies for each bootstrap iteration
Pi.boot <- lapply(res.boot, get_Pi)

### Function returning a list of estimated true trio frequencies (each entry corresponding to a separate MAF interval)
get_Pi_c_list <- function(Pi_parents_df, gamma_list) {
  Pi_parents_list <- split(Pi_parents_df, seq(nrow(Pi_parents_df)))
  Pi_c_list <- map2(Pi_parents_list, gamma_list, get_Pi_c)
  return(Pi_c_list)
}

### Obtain list of estimated true trio frequencies for each bootstrap iteration
Pic.boot <- map2(Pi.boot, gamma.boot, get_Pi_c_list)

### Function returning a list of error rates (each entry corresponding to a separate MAF interval)
get_delta_list <- function(pnt.boot, Pic_list) {
  pnt_list <- split(pnt.boot, seq(nrow(pnt.boot)))
  res <- map2(pnt_list, Pic_list, get_error_rate)
  return(res)
}

### Obtain list of error rates for each bootstrap iteration
delta.boot <- map2(pnt.boot, Pic.boot, get_delta_list)

#####################################################################

### Load in UK Biobank data (this is only used to obtain the number of trio genotypes from each interval used in the weighted average)
ukbio <- read.table(input.file, header = TRUE)
G.df <- ukbio[,c(1,3,4:21)] %>% 
  unlist() %>% 
  matrix(nrow = 20, byrow = TRUE) %>% 
  t() %>%
  as.data.frame() %>%
  filter(V2 != "SUM" & V2 != "EXP" & V2 != "OBS") %>%
  arrange(V1)

### Load in resampled trios for each bootstrap iteration
resamp.1 <- readRDS("/projects/browning/brwnlab/masakin/deletions/Fall_2023/real_data_results/resamp.trios.13.1.rds")
resamp.2 <- readRDS("/projects/browning/brwnlab/masakin/deletions/Fall_2023/real_data_results/resamp.trios.13.2.rds")
resamp.3 <- readRDS("/projects/browning/brwnlab/masakin/deletions/Fall_2023/real_data_results/resamp.trios.13.3.rds")
resamp.4 <- readRDS("/projects/browning/brwnlab/masakin/deletions/Fall_2023/real_data_results/resamp.trios.13.4.rds")
resamp.5 <- readRDS("/projects/browning/brwnlab/masakin/deletions/Fall_2023/real_data_results/resamp.trios.13.5.rds")
resamp.6 <- readRDS("/projects/browning/brwnlab/masakin/deletions/Fall_2023/real_data_results/resamp.trios.13.6.rds")
resamp.7 <- readRDS("/projects/browning/brwnlab/masakin/deletions/Fall_2023/real_data_results/resamp.trios.13.7.rds")
resamp.8 <- readRDS("/projects/browning/brwnlab/masakin/deletions/Fall_2023/real_data_results/resamp.trios.13.8.rds")

### Remove last four bootstrap iterations like before
resamp <- c(resamp.1,
            resamp.2,
            resamp.3,
            resamp.4,
            resamp.5,
            resamp.6,
            resamp.7,
            resamp.8)[-c(101,102,103,104)]

### Returns the average error rate and deletion rate, when the vector of error and deletion rates are from a specific bootstrap iteration with its own set of resampled trios
get_weighted_average <- function(delta_list, gamma_list, trios) {
  ### The number of intervals is two smaller than what is in the dataset (because we're not including (0, 0] and (0, 0.001])
  n_intervals <- G.df %>% filter(V2 == trios[1]) %>% nrow() 
  n_intervals <- n_intervals - 2
  
  ### We accumulate all of the trio genotypes from each trio id
  total <- rep(0, n_intervals*18) %>% matrix(nrow = n_intervals) %>% as.data.frame()
  for (tr in trios) {
    G.filt <- G.df %>% filter(V2 == tr) %>% arrange(V1)
    values <- G.filt[,3:20] %>% 
      mutate_all(function(x) as.numeric(as.character(x))) %>%
      apply(2, undo_cumsum) %>% 
      as.data.frame() %>%
      mutate_all(function(x) as.numeric(as.character(x)))
    values <- values[-c(1,2),]
    total <- add(list(total, values))
  }
  
  ### G.list is split according to the MAF intervals
  G.list <- split(total, seq(nrow(total)))
  
  w <- lapply(G.list, sum) %>% unlist()
  average_delta <- weighted.mean(delta_list %>% unlist(), w)
  average_gamma <- weighted.mean(gamma_list %>% unlist(), w)
  return(c(average_delta, average_gamma))
}

res <- pmap(list(delta.boot, gamma.boot, resamp), get_weighted_average)
saveRDS(res, "weighted.averages.real.data.boot.rds")
