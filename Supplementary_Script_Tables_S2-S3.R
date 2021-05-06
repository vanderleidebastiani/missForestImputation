## Using phylogenetic information to impute missing functional trait values in ecological databases
# Code for results of the tables S2 and S3

## Load packages
require(ape)
require(geiger)
require(picante)
require(phytools)

## Load functions
source("Supplementary_functions.R")

## Set all parameters
# In this example trees are set with 100 species and dataset with 3, 5 and 10 traits (table S3)
PAR1 <- expand.grid(birth = c(0.1), # Speciation rate in the phylogenetic tree
				   death = c(0), # Extinction rate in the phylogenetic tree
				   n.spp = c(100), # Number of species in the phylogenetic tree and in the dataset
				   lambda = c(0.2, 0.6, 1, 1.2), # Phylogenetic signal, Pagel's lambda and Grafen's rho
				   dependent = TRUE, # Logical argument to specify if set correlation between traits 
				   cor.traits = c(0.1, 0.25, 0.5, 0.75, 0.9), # Correlation between traits
				   n.traits = c(3, 5, 10)) # Number of species traits
PAR2 <- expand.grid(birth = c(0.1),
				   death = c(0),
				   n.spp = c(100),
				   lambda = c(0.2, 0.6, 1, 1.2),
				   dependent = FALSE,
				   cor.traits = c(0),
				   n.traits = c(3, 5, 10))

# Merge parameters
PAR <- rbind(PAR1, PAR2)
PAR

## Run the simulation 

# Set main parameters for the simulation
n.rep <- 1000 # Number of repetitions

RESULTS.K.list <- vector("list", nrow(PAR))
for(i in 1:nrow(PAR)){
	print(i)
	parameters <- PAR[i,]
	RESULTS.K.list[[i]] <- RUN.K(parameters, n.rep = n.rep)
}

# Organize the results of K statistic
RESULTS.k <- matrix(NA, nrow(PAR), ncol(PAR)+10)
for(i in 1:nrow(PAR)){
	parameters <- PAR[i,]
	# Organize the results
	results <- RESULTS.K.list[[i]]
	# Summary statistic of results
	results.mean <- apply(results, 2, mean)
	names(results.mean) <- paste0("k_", names(results.mean))
	RESULTS.k[i, ] <- unlist(c(parameters, results.mean, rep(NA, 10-parameters$n.traits)))
	# Set names to RESULTS.k object
	if(i==1){
		colnames(RESULTS.k) <- c(names(parameters), paste0("k_t.", seq_len(10)))
	}
}
## Final results of K statistic
RESULTS.k <- as.data.frame(RESULTS.k)
RESULTS.k

