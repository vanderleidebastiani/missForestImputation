## Using phylogenetic information to impute missing functional trait values in ecological databases
# Code for results of the tables S10 and S11

## Load packages
require(ape)
require(geiger)
require(PVR)
require(phytools)
require(missForest)
require(parallel)

## Load functions
source("Supplementary_functions.R")

## Set all parameters
# In this example trees are set with 100 species and dataset with 3, 5 and 10 traits (table S11)
PAR1 <- expand.grid(birth = c(0.1), # Speciation rate in the phylogenetic tree
					death = c(0), # Extinction rate in the phylogenetic tree
					n.spp = c(100), # Number of species in the phylogenetic tree and in the dataset
					lambda = c(0.2, 0.6, 1, 1.2), # Phylogenetic signal, Pagel lambda and Grafen rho
					n.traits = c(3, 5,10), # Number of species traits
					dependent = TRUE, # Logical to specify if set correlation between traits
					cor.traits = c(0.1, 0.25, 0.5, 0.75, 0.9), # Correlation between traits
					noNA = c(0.05, 0.10, 0.20, 0.30, 0.40, 0.50)) # Proportion of entries cell in trait dataset
PAR2 <- expand.grid(birth = c(0.1), 
					death = c(0), 
					n.spp = c(100), 
					lambda = c(0.2, 0.6, 1, 1.2), 
					n.traits = c(3, 5,10), 
					dependent = FALSE,
					cor.traits = c(0),
					noNA = c(0.05, 0.10, 0.20, 0.30, 0.40, 0.50))
# Merge parameters
PAR <- rbind(PAR1, PAR2)
PAR

# Set main parameters to simulation
n.rep <- 1000 # Number of repetition
parallel <- NULL # Number of parallel processes, it depends on the computer that running the analysis

for(i in 1:nrow(PAR)){
	print(i)
	parameters <- PAR[i,]
	# Run missForest imputation to get the NRMSE
	results.temp <- RUN.random(parameters, n.rep = n.rep, parallel = parallel)
	# Organize the results
	results <- sapply(results.temp, function(x) x$Error, simplify = FALSE)
	results <- do.call(rbind, results)
	rownames(results) <- paste0("list.", seq_len(n.rep), ".nrmse")
	# Base file name
	base <- paste(names(parameters), parameters, sep = "_", collapse = "_")
	# Write the .csv file
	write.csv(results, paste0("./Results/RES_", sprintf("%.4d", i), "_random_nrmse_", base, ".csv"))
}


## Organize the results

# Small function to get summary statistics
f2 <- function(x){
	r.mean.temp <- apply(x, 2, mean)
	r.sd.temp <- apply(x, 2, sd)
	res.temp <- c(r.mean.temp[1], r.sd.temp[1]) # Mean and standard deviation to NRMSE with Random imputation
	names(res.temp) <- c("mean.nrmse", 
						 "sd.nrmse")
	return(res.temp)
}

# Organize the results of NRMSE
RESULTS.nrmse.random <- matrix(NA, nrow(PAR), ncol(PAR)+2)
for(i in 1:nrow(PAR)){
	print(i)
	parameters <- PAR[i,]
	# Base file name
	base <- paste(names(parameters), parameters, sep = "_", collapse = "_")
	# Read .csv with the results
	res.nrmse <- read.csv(paste0("./Results/RES_", sprintf("%.4d", i), "_random_nrmse_", base, ".csv"), row.names = 1)
	# Summary statistic of results
	res.temp.nrmse.summary <- f2(res.nrmse)
	RESULTS.nrmse.random[i, ] <- unlist(c(parameters, res.temp.nrmse.summary))
	# Set names to RESULTS.nrmse object
	if(i==1){
		colnames(RESULTS.nrmse.random) <- names(unlist(c(parameters, res.temp.nrmse.summary)))
	}
}

## Final results to NRMSE with Random imputation
RESULTS.nrmse.random <- as.data.frame(RESULTS.nrmse.random)
RESULTS.nrmse.random

