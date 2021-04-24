## Using phylogenetic information to impute missing functional trait values in ecological databases
# Code for results of the table S1

## Load packages
require(ape)
require(geiger)
require(picante)
require(phytools)

## Load functions
source("Supplementary_functions.R")

## Set all parameters
PAR <- expand.grid(birth = c(0.1), # Speciation rate in the phylogenetic tree
				   death = c(0), # Extinction rate in the phylogenetic tree
				   n.spp = c(50, 100), # Number of species in the phylogenetic tree and in the dataset
				   lambda = c(0.2, 0.6, 1, 1.2), # Phylogenetic signal, Pagel lambda and Grafen rho
				   dependent = FALSE, # Logical to specify if set correlation between traits 
				   cor.traits = c(0), # Correlation between traits
				   n.traits = c(1)) # Number of species traits
PAR


## Run the simulation 

# Set main parameters to simulation
n.rep <- 1000 # Number of repetition

RESULTS.K.list <- vector("list", nrow(PAR))
for(i in 1:nrow(PAR)){
	print(i)
	parameters <- PAR[i,]
	RESULTS.K.list[[i]] <- RUN.K(parameters, n.rep = n.rep)
}

# Organize the results of K statistic
RESULTS.k <- matrix(NA, nrow(PAR), ncol(PAR)+2)
for(i in 1:nrow(PAR)){
	print(i)
	parameters <- PAR[i,]
	# Organize the results
	results <- t(RESULTS.K.list[[i]])
	colnames(results) <- paste0("k")
	# Summary statistic of results
	results.mean.sd <- c(apply(results, 2, mean), 
						 apply(results, 2, sd))
	names(results.mean.sd) <- c("mean", "sd")
	RESULTS.k[i, ] <- unlist(c(parameters, results.mean.sd))
	# Set names to RESULTS.k object
	if(i==1){
		colnames(RESULTS.k) <- names(unlist(c(parameters, results.mean.sd)))
	}
}
## Final results of K statistic
RESULTS.k <- as.data.frame(RESULTS.k)
RESULTS.k

