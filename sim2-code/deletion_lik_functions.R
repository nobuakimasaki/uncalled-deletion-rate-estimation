##### Functions used to obtain genotype frequencies and likelihood of model

### Returns estimate of ancestral parental genotype frequencies using Mendelian-consistent trio genotype counts
get_Pi_parents_consistent <- function(G) {
  G000 <- G[1]
  G001 <- G[2]
  G002 <- G[3]
  G010 <- G[4]
  G011 <- G[5]
  G012 <- G[6] 
  G020 <- G[7]
  G021 <- G[8]
  G022 <- G[9]
  G110 <- G[10]
  G111 <- G[11]  
  G112 <- G[12]
  G120 <- G[13]
  G121 <- G[14]
  G122 <- G[15]
  G220 <- G[16]
  G221 <- G[17]
  G222 <- G[18]
  n <- sum(G000, G010, G011, G021, G110, G111, G112, G121, G122, G222)
  Pi00 <- (G000)/n
  Pi01 <- (G010+G011)/n
  Pi02 <- (G021)/n
  Pi11 <- (G110+G111+G112)/n
  Pi12 <- (G121+G122)/n
  Pi22 <- (G222)/n
  return(c(Pi00, Pi01, Pi02, Pi11, Pi12, Pi22))
}

### Returns estimate of true trio genotype frequencies using estimate of ancestral parental genotype frequencies and the deletion rate 
get_Pi_c <- function(Pi_parents, gamma) {
  r <- 1-2*gamma
  c <- gamma
  Pi00 <- Pi_parents[1] 
  Pi01 <- Pi_parents[2] 
  Pi02 <- Pi_parents[3]
  Pi11 <- Pi_parents[4]
  Pi12 <- Pi_parents[5]
  Pi22 <- Pi_parents[6]
  
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

### Term inside summation in Equation 1
Pi_theta_term <- function(i,j,k,theta,Pi,x,y,z) {
  if (x != y) {
    term <- Pi[i,j,k]*(theta[i,x]*theta[j,y] + theta[i,y]*theta[j,x])*theta[k,z]
  }
  else {
    term <- Pi[i,j,k]*theta[i,x]*theta[j,y]*theta[k,z]
  }
  return(term)
}

### Equation 1
lik <- function(x,y,z,theta,Pi) {
  l <- data.frame(l=1:5)
  m <- data.frame(m=1:5)
  n <- data.frame(n=1:5)
  d_ <- merge(l,m,all=TRUE) %>% merge(n,all=TRUE)
  d_ <- d_ %>% filter(m >= l) %>% arrange(l,m,n)
  l <- d_$l
  m <- d_$m
  n <- d_$n
  
  summed <- pmap(list(l,m,n), Pi_theta_term, theta=theta, Pi=Pi, x=x, y=y, z=z) %>%
    unlist() %>%
    sum()
  return(summed)
}

### Equation 2
neg_sum_log_lik <- function(theta01,theta02,theta10,theta12,theta20,theta21,theta31,theta32,theta40,theta41,gamma,Pi_parents,G) {
  theta00 <- 1-theta01-theta02
  theta11 <- 1-theta10-theta12
  theta22 <- 1-theta20-theta21
  theta30 <- 1-theta31-theta32
  theta42 <- 1-theta40-theta41
  
  pen <- 0

  for (param in c(theta00, theta11, theta22, theta30, theta42)) {
    if (param < 0) {
      pen <- pen + param^2
    }
  }
  if (pen > 0) {return(1e+13 + pen)}
  
  if (sum(is.na(Pi_parents)) > 0 | sum(Pi_parents < 0) > 0) {return(1e+15)}
  
  theta <- c(theta00, theta01, theta02,
             theta10, theta11, theta12,
             theta20, theta21, theta22,
             theta30, theta31, theta32,
             theta40, theta41, theta42)
  theta <- matrix(theta,
                  nrow = 5,
                  byrow = TRUE)
  Pic <- get_Pi_c(Pi_parents, gamma)
  
  x <- data.frame(x=1:3)
  y <- data.frame(y=1:3)
  z <- data.frame(z=1:3)
  d_ <- merge(x,y,all=TRUE) %>% merge(z,all=TRUE)
  d_ <- d_ %>% dplyr::filter(y >= x) %>% arrange(x,y,z)
  x <- d_$x
  y <- d_$y
  z <- d_$z
  
  log_lik <- pmap(list(x,y,z),lik,theta=theta,Pi=Pic) %>% unlist() %>% log()
  sum_log_lik <- sum(G * log_lik)
  return(-sum_log_lik)
}

# get_Pi_parents <- function(G) {
#   G000 <- G[1]
#   G001 <- G[2]
#   G002 <- G[3]
#   G010 <- G[4]
#   G011 <- G[5]
#   G012 <- G[6] 
#   G020 <- G[7]
#   G021 <- G[8]
#   G022 <- G[9]
#   G110 <- G[10]
#   G111 <- G[11]  
#   G112 <- G[12]
#   G120 <- G[13]
#   G121 <- G[14]
#   G122 <- G[15]
#   G220 <- G[16]
#   G221 <- G[17]
#   G222 <- G[18]
#   n <- sum(G)
#   Pi00 = (G000 + G001 + G002) / n
#   Pi01 = (G010 + G011 + G012) / n
#   Pi02 = (G020 + G021 + G022) / n
#   Pi11 = (G110 + G111 + G112) / n
#   Pi12 = (G120 + G121 + G122) / n
#   Pi22 = (G220 + G221 + G222) / n
#   return(c(Pi00, Pi01, Pi02, Pi11, Pi12, Pi22))
# }
# get_Pi <- function(Pi_parents) {
#   Pi00 <- Pi_parents[1]
#   Pi01 <- Pi_parents[2]
#   Pi02 <- Pi_parents[3]
#   Pi11 <- Pi_parents[4]
#   Pi12 <- Pi_parents[5]
#   Pi22 <- Pi_parents[6]
#   
#   Pi000 <- Pi00
#   Pi010 <- Pi01/2
#   Pi011 <- Pi01/2
#   Pi021 <- Pi02
#   
#   Pi110 <- Pi11/4
#   Pi111 <- Pi11/2
#   Pi112 <- Pi11/4
#   Pi121 <- Pi12/2
#   Pi122 <- Pi12/2
#   
#   Pi222 <- Pi22
#   
#   Pi0 <- matrix(c(Pi000, 0, 0, Pi010, Pi011, 0, 0, Pi021, 0), nrow = 3, ncol = 3, byrow = TRUE)
#   Pi1 <- matrix(c(0, 0, 0, Pi110, Pi111, Pi112, 0, Pi121, Pi122), nrow = 3, ncol = 3, byrow = TRUE)
#   Pi2 <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, Pi222), nrow = 3, ncol = 3, byrow = TRUE)
#   
#   Pi <- array(c(Pi0, Pi1, Pi2), dim = c(3,3,3))
#   ## Pi^{ijk} is Pi[j,k,i]
#   
#   Pi <- aperm(Pi, c(3,1,2))
#   ## Pi^{ijk} is Pi[i,j,k]
#   return(Pi)
# }

