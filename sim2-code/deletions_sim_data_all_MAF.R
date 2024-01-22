##### Functions used to simulate trio genotypes for 100 MAF intervals

set.seed(25)
library(dplyr)
library(purrr)
library(parallel)
source("deletion_sim_functions_all_MAF.R")
source("deletion_lik_functions.R")
n.cores <- detectCores()

### Simulates trio genotypes and returns point estimates and bootstrapped estimates
### n_boot: number of bootstrap iterations
sim_data_all_MAF <- function(n_trios, n_pos_vec, p_mat, theta, gamma) {
  n_reps <- 10
  n_pos_vec <- n_pos_vec/10
  ### Simulate trio genotype counts
  res <- pmap(list(p_mat[,1], p_mat[,2], n_pos_vec), gen_dataset,
             n_trios = n_trios, 
             theta = theta, 
             gamma = gamma,
             n_reps)
  return(res)
}

p0 <- c(0, 0.001)
p1 <- c(0.001, 0.005)
p_matrix <- matrix(c(p0, p1), nrow = 2, byrow = TRUE)
p_start <- seq(0.005, 0.495, 0.005)
p_matrix_add <- cbind(p_start, p_start + 0.005)
p_matrix <- rbind(p_matrix, p_matrix_add)

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

### trio genotype counts we want to mimic
n_trio_genotypes_vec <- c(
  32960735766, 401791703, 91112338, 44170814, 30130478, 22841778,
  18472996, 15706022, 13822779, 12167986, 11034631, 10091069,
  9578837, 8930988, 8308059, 7889474, 7361688, 7171047,
  6899507, 6370733, 6412560, 6144866, 6002204, 5844377,
  5763633, 5602302, 5517174, 5263140, 5170935, 5070166,
  4958062, 4825304, 4744784, 4664949, 4597961, 4505360,
  4456003, 4528612, 4279702, 4237813, 4191658, 4161193,
  4103316, 4019117, 4160791, 3944178, 3852143, 3880736,
  3747152, 3786694, 3755833, 3689185, 3663670, 3617182,
  3667778, 3688273, 3511297, 3407225, 3434725, 3492256,
  3443026, 3303085, 3431335, 3353170, 3373477, 3333397,
  3418652, 3305150, 3245374, 3197937, 3263587, 3236900,
  3213430, 3200723, 3214038, 3241592, 3182956, 3185683,
  3157201, 3091805, 3070250, 3092941, 3030813, 3090433,
  3059397, 3056895, 3013304, 3030130, 3055726, 3049241,
  2995689, 3060875, 2950950, 3031150, 3013575, 3047227,
  3079978, 3063788, 2967030, 2954832, 2991574
)

### round to largest digit and divide by 100
n_pos_vec_rounded <- round(n_trio_genotypes_vec, -floor(log10(n_trio_genotypes_vec)))/100

# res <- sim_data_all_MAF(100, n_pos_vec_rounded[2:100], p_matrix[2:100,], theta, gamma)
# debugonce(sim_data_all_MAF)
# debugonce(gen_dataset)
# debugonce(gen_data)
# debugonce(get_G_sim2)

### generate trio genotype counts and put into dataframe
res <- sim_data_all_MAF(100, n_pos_vec_rounded[2:101], p_matrix[2:101,], theta, gamma)

trio_counts <- matrix(nrow = 0, ncol = 4) %>% as.data.frame()
colnames(trio_counts) <- c("par_1", "par_2", "offsp", "count")
for (i in 1:length(res)) {
  temp <- res[[i]][[1]] %>% as.data.frame()
  temp$trio_id <- 1:100
  temp$MAF_start <- p_matrix[i+1,1]
  trio_counts <- rbind(trio_counts, temp)
}

### the second argument contains the pre-deletion parental genotype counts
c_parents_anc <- res %>% lapply(function(x) {x[[2]]})
### the third argument contains the true trio genotype counts
c_trios_true <- res %>% lapply(function(x) {x[[3]]})

saveRDS(trio_counts, "sim_obs_trio_cnts_all_MAF.rds")
saveRDS(c_parents_anc, "sim_predel_parental_cnts_all_MAF.rds")
saveRDS(c_trios_true, "sim_true_trio_cnts_all_MAF.rds")

