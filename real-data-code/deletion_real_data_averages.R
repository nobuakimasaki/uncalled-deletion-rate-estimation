##### Code to obtain point estimates for the average genotype error rate and deletion rate. Requires the file of point estimates to be in the same directory. 

### Load in packages
library(dplyr)
library(tibble)
library(tidyr)
library(purrr)

args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("Input file must be supplied", call.=FALSE) }

input.file <- args[1]

### Helper functions to add matrices together and undo cumulative sum
add <- function(x) Reduce("+", x)

undo_cumsum <- function(x) {
  c(x[1],diff(x))
}

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

### Function returning a matrix of point estimates, each row representing a separate MAF interval
### i: index specifying model and optimization algorithm used
get_est <- function(pnt, i) {
  est.df <- lapply(pnt, function(x) {x[[i]][[1]]}) %>% unlist() %>% matrix(nrow = length(pnt), byrow = TRUE)  %>% exp() %>% as.data.frame()
  return(est.df)
}

### Function returning a matrix of estimated ancestral parental frequencies, each row representing a separate MAF interval
get_Pi <- function(pnt) {
  est.df <- lapply(pnt, function(x) {x[[1]]}) %>% unlist() %>% matrix(nrow = length(pnt), byrow = TRUE) %>% as.data.frame()
  return(est.df)
}

### Function returning error rate for a single MAF interval
### pnt: vector of point estimates
get_error_rate <- function(pnt, Pic) {
  Pik_vec <- rep(NA, 5)
  # sum together trio frequencies to get frequency of each genotype
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
  
  # if genotype does not contain deleted allele, multiply by probability of miscall
  delta <- 0
  delta_from_deletion <- 0
  for (k in 1:3) {
    delta <- delta + Pik_vec[k]*(1-theta_vec[k])
  }
  # if genotype contains deleted allele, it is miscalled
  for (k in 4:5) {
    delta <- delta + Pik_vec[k]
    delta_from_deletion <- delta_from_deletion + Pik_vec[k]
  }
  
  return(c(delta, delta_from_deletion))
}

### Returns weighted average of deletion rate and genotype error rate
### delta_list: vector of genotype error rates (each entry representing the overall error rate and the error rate due to deletions for a single MAF interval)
get_weighted_average <- function(delta_list, gamma_list, G.list) {
  w <- lapply(G.list, sum) %>% unlist()
  print("weights used: ")
  print(w)
  delta_from_all_list <- delta_list %>% lapply(function(x) {x[1]}) 
  delta_from_deletion_list <- delta_list %>% lapply(function(x) {x[2]}) 
  
  average_delta <- weighted.mean(delta_from_all_list %>% unlist(), w)
  average_delta_from_deletion <- weighted.mean(delta_from_deletion_list %>% unlist(), w)
  
  average_gamma <- weighted.mean(gamma_list %>% unlist(), w)
  
  return(c(average_delta, average_delta_from_deletion, average_gamma))
}

### Load in matrix of point estimates
res <- readRDS("/projects/browning/brwnlab/masakin/deletions/Fall_2023/real_data_results/est.res.pt.rds")
#res <- readRDS("result/est.res.pt.rds")

### Organize point estimates and remove first row (MAF interval: (0, 0.001])
pnt <- get_est(res, 6)
colnames(pnt) <- c("t01", "t02", "t10", "t12", "t20", "t21", "g")
pnt.list <- split(pnt, 1:nrow(pnt))
pnt.list <- pnt.list[-1]
### Obtain estimated deletion rate for each MAF interval
gamma_list <- pnt$g %>% as.list()
gamma_list <- gamma_list[-1]

### Run previously defined functions to get estimated error rate for each MAF interval
Pi_parents <- get_Pi(res)
Pi_parents_list <- split(Pi_parents, 1:nrow(Pi_parents))
Pi_parents_list <- Pi_parents_list[-1]
Pic_list <- map2(Pi_parents_list, gamma_list, get_Pi_c)
delta_list <- map2(pnt.list, Pic_list, get_error_rate)

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
### Get unique trio ids
trios <- unique(G.df$V2)

### The number of intervals is two smaller than what is in this dataset (because we're not including (0, 0] and (0, 0.001])
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
  values <- values[-c(1:2),]
  total <- add(list(total, values))
}
### G.list is split according to the MAF intervals
G.list <- split(total, seq(nrow(total)))
print("UK Biobank counts: ")
print(G.list)

### Get weighted averages and save
res <- get_weighted_average(delta_list, gamma_list, G.list)
saveRDS(res, "weighted.averages.real.data.rds")
