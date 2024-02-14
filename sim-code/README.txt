This directory stores the files used to simulate the data and apply the model in the first simulation study. 

deletion_analysis.rmd is the file used to analyze the results of fitting the model to this data.

deletion_lik_functions.R contains the functions used to calculate the likelihood when fitting the model.

deletion_sim_functions.R contains the functions used to simulate the observed trio genotypes.

Each of the deletions_sim_cluster_*.R files represents the script used to generate the observed trio genotypes from a specific MAF interval and fit the model to the generated data.
