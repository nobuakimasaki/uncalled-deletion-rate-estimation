##### Functions used to simulate trio genotypes for simulation study

### Returns a genotype vector with deletions according to the deletion rate
add_deletions <- function(genotype, gamma) {
  genotype[genotype == "AA"] <- sample(c("AA", "AD"), length(genotype[genotype == "AA"]), replace = TRUE, prob = c(1-2*gamma, 2*gamma))
  genotype[genotype == "AB"] <- sample(c("AB", "AD", "BD"), length(genotype[genotype == "AB"]), replace = TRUE, prob = c(1-2*gamma, gamma, gamma))
  genotype[genotype == "BB"] <- sample(c("BB", "BD"), length(genotype[genotype == "BB"]), replace = TRUE, prob = c(1-2*gamma, 2*gamma))
  return(genotype)
}

### Returns a genotype vector of the offspring based on transmission equilibrium
get_offsp <- function(parent1, parent2) {
  offsp <- rep(NA, length(parent1))
  offsp[parent1 == "AA" & parent2 == "AA"] <- "AA"
  
  offsp[parent1 == "AA" & parent2 == "AB" | parent1 == "AB" & parent2 == "AA"] <- 
    sample(c("AA", "AB"), length(offsp[parent1 == "AA" & parent2 == "AB" | parent1 == "AB" & parent2 == "AA"]), replace = TRUE, prob = c(0.5, 0.5))
  
  offsp[parent1 == "AA" & parent2 == "BB" | parent1 == "BB" & parent2 == "AA"] <- "AB"
  
  offsp[parent1 == "AA" & parent2 == "AD" | parent1 == "AD" & parent2 == "AA"] <- 
    sample(c("AA", "AD"), length(offsp[parent1 == "AA" & parent2 == "AD" | parent1 == "AD" & parent2 == "AA"]), replace = TRUE, prob = c(0.5, 0.5))
  
  offsp[parent1 == "AA" & parent2 == "BD" | parent1 == "BD" & parent2 == "AA"] <- 
    sample(c("AB", "AD"), length(offsp[parent1 == "AA" & parent2 == "BD" | parent1 == "BD" & parent2 == "AA"]), replace = TRUE, prob = c(0.5, 0.5))
  
  offsp[parent1 == "AB" & parent2 == "AB"] <- sample(c("AA", "AB", "BB"), length(offsp[parent1 == "AB" & parent2 == "AB"]), replace = TRUE, prob = c(0.25, 0.5, 0.25))
  
  offsp[parent1 == "AB" & parent2 == "BB" | parent1 == "BB" & parent2 == "AB"] <- 
    sample(c("AB", "BB"), length(offsp[parent1 == "AB" & parent2 == "BB" | parent1 == "BB" & parent2 == "AB"]), replace = TRUE, prob = c(0.5, 0.5))
  
  offsp[parent1 == "AB" & parent2 == "AD" | parent1 == "AD" & parent2 == "AB"] <- 
    sample(c("AA", "AB", "AD", "BD"), length(offsp[parent1 == "AB" & parent2 == "AD" | parent1 == "AD" & parent2 == "AB"]), replace = TRUE, prob = c(0.25, 0.25, 0.25, 0.25))
  
  offsp[parent1 == "AB" & parent2 == "BD" | parent1 == "BD" & parent2 == "AB"] <- 
    sample(c("AB", "BB", "AD", "BD"), length(offsp[parent1 == "AB" & parent2 == "BD" | parent1 == "BD" & parent2 == "AB"]), replace = TRUE, prob = c(0.25, 0.25, 0.25, 0.25))
  
  offsp[parent1 == "BB" & parent2 == "BB"] <- "BB"
  
  offsp[parent1 == "BB" & parent2 == "AD" | parent1 == "AD" & parent2 == "BB"] <- 
    sample(c("AB", "BD"), length(offsp[parent1 == "BB" & parent2 == "AD" | parent1 == "AD" & parent2 == "BB"]), replace = TRUE, prob = c(0.5, 0.5))
  
  offsp[parent1 == "BB" & parent2 == "BD" | parent1 == "BD" & parent2 == "BB"] <- 
    sample(c("BB", "BD"), length(offsp[parent1 == "BB" & parent2 == "BD" | parent1 == "BD" & parent2 == "BB"]), replace = TRUE, prob = c(0.5, 0.5))
  
  offsp[parent1 == "AD" & parent2 == "AD"] <- sample(c("AA", "AD", "DD"), length(offsp[parent1 == "AD" & parent2 == "AD"]), replace = TRUE, prob = c(0.25, 0.5, 0.25))
  
  offsp[parent1 == "AD" & parent2 == "BD" | parent1 == "BD" & parent2 == "AD"] <- 
    sample(c("AB", "AD", "BD", "DD"), length(offsp[parent1 == "AD" & parent2 == "BD" | parent1 == "BD" & parent2 == "AD"]), replace = TRUE, prob = c(0.25, 0.25, 0.25, 0.25))
  
  offsp[parent1 == "BD" & parent2 == "BD"] <- sample(c("BB", "BD", "DD"), length(offsp[parent1 == "BD" & parent2 == "BD"]), replace = TRUE, prob = c(0.25, 0.5, 0.25))
  return(offsp)
}

### Returns a genotype vector with simulated genotype errors 
sim_err <- function(genotype, theta) {
  genotype[genotype == "AA"] <- sample(c("AA", "AB", "BB"),
                                       length(genotype[genotype == "AA"]),
                                       replace = TRUE,
                                       prob = c(theta[1],
                                                theta[2],
                                                theta[3]))
  genotype[genotype == "AB"] <- sample(c("AA", "AB", "BB"),
                                       length(genotype[genotype == "AB"]),
                                       replace = TRUE,
                                       prob = c(theta[4],
                                                theta[5],
                                                theta[6]))
  genotype[genotype == "BB"] <- sample(c("AA", "AB", "BB"),
                                       length(genotype[genotype == "BB"]),
                                       replace = TRUE,
                                       prob = c(theta[7],
                                                theta[8],
                                                theta[9]))
  genotype[genotype == "AD"] <- sample(c("AA", "AB", "BB"),
                                       length(genotype[genotype == "AD"]),
                                       replace = TRUE,
                                       prob = c(theta[10],
                                                theta[11],
                                                theta[12]))
  genotype[genotype == "BD"] <- sample(c("AA", "AB", "BB"),
                                       length(genotype[genotype == "BD"]),
                                       replace = TRUE,
                                       prob = c(theta[13],
                                                theta[14],
                                                theta[15]))
  return(genotype)
}

### Counts the number of each trio genotype
get_G_sim <- function(parent1, parent2, offsp) {
  G000 <- sum(parent1 == "AA" & parent2 == "AA" & offsp == "AA")
  G001 <- sum(parent1 == "AA" & parent2 == "AA" & offsp == "AB")
  G002 <- sum(parent1 == "AA" & parent2 == "AA" & offsp == "BB")
  
  G010 <- sum(parent1 == "AA" & parent2 == "AB" & offsp == "AA" | parent1 == "AB" & parent2 == "AA" & offsp == "AA")
  G011 <- sum(parent1 == "AA" & parent2 == "AB" & offsp == "AB" | parent1 == "AB" & parent2 == "AA" & offsp == "AB")
  G012 <- sum(parent1 == "AA" & parent2 == "AB" & offsp == "BB" | parent1 == "AB" & parent2 == "AA" & offsp == "BB")
  
  G020 <- sum(parent1 == "AA" & parent2 == "BB" & offsp == "AA" | parent1 == "BB" & parent2 == "AA" & offsp == "AA")
  G021 <- sum(parent1 == "AA" & parent2 == "BB" & offsp == "AB" | parent1 == "BB" & parent2 == "AA" & offsp == "AB")
  G022 <- sum(parent1 == "AA" & parent2 == "BB" & offsp == "BB" | parent1 == "BB" & parent2 == "AA" & offsp == "BB")
  
  G110 <- sum(parent1 == "AB" & parent2 == "AB" & offsp == "AA")
  G111 <- sum(parent1 == "AB" & parent2 == "AB" & offsp == "AB")
  G112 <- sum(parent1 == "AB" & parent2 == "AB" & offsp == "BB")
  
  G120 <- sum(parent1 == "AB" & parent2 == "BB" & offsp == "AA" | parent1 == "BB" & parent2 == "AB" & offsp == "AA")
  G121 <- sum(parent1 == "AB" & parent2 == "BB" & offsp == "AB" | parent1 == "BB" & parent2 == "AB" & offsp == "AB")
  G122 <- sum(parent1 == "AB" & parent2 == "BB" & offsp == "BB" | parent1 == "BB" & parent2 == "AB" & offsp == "BB")
  
  G220 <- sum(parent1 == "BB" & parent2 == "BB" & offsp == "AA")
  G221 <- sum(parent1 == "BB" & parent2 == "BB" & offsp == "AB")
  G222 <- sum(parent1 == "BB" & parent2 == "BB" & offsp == "BB")
  
  G <- c(G000, G001, G002, G010, G011, G012, G020, G021, G022,
         G110, G111, G112, G120, G121, G122,
         G220, G221, G222)
  return(G)
}

### Calculates parental genotype frequencies
get_Pi_parents_sim <- function(parent1, parent2) {
  n <- length(parent1)
  Pi00 <- sum(parent1 == "AA" & parent2 == "AA")/n
  Pi01 <- sum(parent1 == "AA" & parent2 == "AB" | parent1 == "AB" & parent2 == "AA")/n
  Pi02 <- sum(parent1 == "AA" & parent2 == "BB" | parent1 == "BB" & parent2 == "AA")/n
  
  Pi11 <- sum(parent1 == "AB" & parent2 == "AB")/n
  Pi12 <- sum(parent1 == "AB" & parent2 == "BB" | parent1 == "BB" & parent2 == "AB")/n
  
  Pi22 <- sum(parent1 == "BB" & parent2 == "BB")/n
  return(c(Pi00, Pi01, Pi02, Pi11, Pi12, Pi22))
}

### Calculates parental genotype counts
get_Pi_parents_count_sim <- function(parent1, parent2) {
  c00 <- sum(parent1 == "AA" & parent2 == "AA")
  c01 <- sum(parent1 == "AA" & parent2 == "AB" | parent1 == "AB" & parent2 == "AA")
  c02 <- sum(parent1 == "AA" & parent2 == "BB" | parent1 == "BB" & parent2 == "AA")
  
  c11 <- sum(parent1 == "AB" & parent2 == "AB")
  c12 <- sum(parent1 == "AB" & parent2 == "BB" | parent1 == "BB" & parent2 == "AB")
  
  c22 <- sum(parent1 == "BB" & parent2 == "BB")
  return(c(c00, c01, c02, c11, c12, c22))
}

### Simulate trio genotypes
### n_trios: number of trios
### n_pos: number of positions to draw genotypes
### p: start of MAF interval used for markers (p, p + 0.005)
gen_data <- function(x, n_trios, n_pos, p, theta, gamma) {
  start <- Sys.time()
  n <- n_trios*n_pos
  
  ### First generate a MAF for each marker
  p.vec <- runif(n_pos, min = p, max = p + 0.005)
  
  ### For each position, generate parents' genotypes based on MAF
  parent1 <- lapply(p.vec, function(p) {sample(c("AA","AB","BB"), size = n_trios, replace = TRUE, prob = c((1-p)^2,2*p*(1-p),p^2))}) %>% unlist()
  parent2 <- lapply(p.vec, function(p) {sample(c("AA","AB","BB"), size = n_trios, replace = TRUE, prob = c((1-p)^2,2*p*(1-p),p^2))}) %>% unlist()
  ### Count number of each ancestral parental genotype
  c_parents_anc <- get_Pi_parents_count_sim(parent1, parent2)
  
  ### At this point, genotypes are clustered by position (i.e. individual1 at pos1, individual2 at pos1, ...)
  
  ### Add deletions and get offspring for each position and parent pair
  parent1.mut <- add_deletions(parent1, gamma = gamma) %>% unlist()
  parent2.mut <- add_deletions(parent2, gamma = gamma) %>% unlist()
  offsp <- get_offsp(parent1.mut, parent2.mut)
  
  ### Add genotype errors and split the trio genotypes by trio instead of position 
  ### (each row represents a position, each column represents a trio)
  parent1.mut <- sim_err(parent1.mut, theta) %>% matrix(byrow = TRUE, ncol = n_trios)
  parent2.mut <- sim_err(parent2.mut, theta) %>% matrix(byrow = TRUE, ncol = n_trios)
  offsp <- sim_err(offsp, theta) %>% matrix(byrow = TRUE, ncol = n_trios)
  
  ### Put each trio in separate item in list, and get trio genotype counts for each trio
  parent1.mut.lst <- lapply(seq_len(ncol(parent1.mut)), function(i) parent1.mut[,i])
  parent2.mut.lst <- lapply(seq_len(ncol(parent2.mut)), function(i) parent2.mut[,i])
  offsp.lst <- lapply(seq_len(ncol(offsp)), function(i) offsp[,i])
  G_list <- pmap(list(parent1.mut.lst,
                      parent2.mut.lst,
                      offsp.lst), get_G_sim) %>% 
    unlist() %>%
    matrix(nrow = n_trios, byrow = TRUE)
  print("time spent data generation: ")
  print( Sys.time() - start )
  print("c parents anc: ")
  print(c_parents_anc)
  return(list(G_list, c_parents_anc))
}

### Function to repeat above n_reps times (due to memory constraints). When using this function, the number of positions is effectively n_pos*n_reps.
gen_dataset <- function(n_trios, n_pos, n_reps, p, theta, gamma) {
  G_list_list <- lapply(1:n_reps, gen_data, n_trios = n_trios, n_pos = n_pos, p = p, theta = theta, gamma = gamma) 
  G_list <- Reduce(`+`, lapply(G_list_list, function(x) x[[1]]))
  c_parents_anc <- lapply(G_list_list, function(x) x[[2]]) %>% unlist() %>% 
    matrix(byrow = TRUE, ncol = 6)
  c_parents_anc_sum <- colSums(c_parents_anc)
  print("c parents anc: ")
  print(c_parents_anc)
  ### We calculate the ancestral parental frequencies, but we don't use them when fitting the model (mainly to see if model performs better if the true frequencies are used)
  Pi_parents_anc <- c_parents_anc_sum/sum(c_parents_anc_sum)
  return(list(G_list, Pi_parents_anc))
}
