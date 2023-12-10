##### We may want to edit this code to remove the first two MAF intervals
##### Code to obtain number of missing genotype trios.

### Load in packages
library(dplyr)
library(purrr)
library(parallel)

### Helper functions to add matrices together and undo cumulative sum
add <- function(x) Reduce("+", x)

undo_cumsum <- function(x) {
  c(x[1],diff(x))
}

### Load data
input.file <- "/projects/browning/ukbio/masakin/chr1-22.filt45.all.nonmaj1.wb.triogt.cnts"
ukbio <- read.table(input.file, header = TRUE)
G.df <- ukbio[,c(1,3,4:22)] %>% 
  unlist() %>% 
  matrix(nrow = 21, byrow = TRUE) %>% 
  t() %>%
  as.data.frame() %>%
  filter(V2 != "SUM" & V2 != "EXP" & V2 != "OBS") %>%
  arrange(V1)

### Get unique trio ids
trios <- unique(G.df$V2)

n_intervals <- G.df %>% filter(V2 == trios[1]) %>% nrow() 
n_intervals <- n_intervals - 1

### We accumulate all of the trio genotypes from each trio id
total <- rep(0, n_intervals*19) %>% matrix(nrow = n_intervals) %>% as.data.frame()
for (tr in trios) {
  G.filt <- G.df %>% filter(V2 == tr) %>% arrange(V1)
  values <- G.filt[,3:21] %>% 
    mutate_all(function(x) as.numeric(as.character(x))) %>%
    apply(2, undo_cumsum) %>% 
    as.data.frame() %>%
    mutate_all(function(x) as.numeric(as.character(x)))
  values <- values[-1,]
  total <- add(list(total, values))
}

### G.list is split according to the MAF intervals
G.list <- split(total, seq(nrow(total)))

print(G.list)
#w <- lapply(G.list, function(x) {sum(x[-length(x)])}) %>% unlist()
#print(w)
#print(sum(w))
print(colSums(total))

print(rowSums(total[,-19]))
print(sum(rowSums(total[,-19])))