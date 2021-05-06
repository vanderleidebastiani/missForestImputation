## Using phylogenetic information to impute missing functional trait values in ecological databases

# Supplementary functions

# Function for data simulation and to get the cumulative proportion of eigenvalue explanation in PVR
# parameters = A data.frame with only one row and three columns, the names of columns are all parameters setting to run the simulation
get.PVR <- function(parameters){
	# Simulate the phylogenetic tree
	tree <- geiger::drop.extinct(geiger::sim.bdtree(b = parameters$birth, d = parameters$death, stop = "taxa", n = parameters$n.spp, extinct = FALSE))
	tree$tip.label <- sprintf("Sp_%.4d", 1:length(tree$tip.label))	
	tree <- ape::makeNodeLabel(tree, method = "number", prefix = "Node")
	# Decomposing phylogenetic distance matrix derived from tree into a set of orthogonal vectors
	x <- PVR::PVRdecomp(tree)
	# Proportion of eigenvalue explanation
	PVR.eigen.prop <- x@Eigen$values/sum(x@Eigen$values)
	# Cumulative proportion of eigenvalue explanation (number species -2)
	Eigen.prop <- cumsum(PVR.eigen.prop)[1:(parameters$n.spp-2)]
	return(Eigen.prop)
}

# Function to repeat the simulation with each parameter combination (PVR)
# parameters = A data.frame with only one row and three columns, the names of columns are all parameters setting to run the simulation
# n.rep = The number of repetition
RUN.PVR <- function(parameters, n.rep){
	f <- function(x,parameters){
		RES <- get.PVR(parameters)
		return(RES)
	}
	RES <- sapply(as.list(1:n.rep), f, parameters = parameters, simplify = FALSE)
	return(RES)
}


# Function to data simulation and get K values for each trait
# parameters = A data.frame with only one row and seven columns, the names of columns are all parameters setting to run the simulation
get.K <- function(parameters){
	# Simulate the phylogenetic tree
	tree <- geiger::drop.extinct(geiger::sim.bdtree(b = parameters$birth, d = parameters$death, stop = "taxa", n = parameters$n.spp, extinct = FALSE))
	tree$tip.label <- sprintf("Sp_%.4d", 1:length(tree$tip.label))	
	tree <- ape::makeNodeLabel(tree, method = "number", prefix = "Node")
	# Transform tree using Pagel's lambda and Grafen's rho methods
	if(parameters$lambda<=1){
		tree.lambda <- geiger::rescale(tree, model = "lambda", parameters$lambda)
	} else{
		tree.lambda <- ape::compute.brlen(tree, method = "Grafen", power = parameters$lambda)
	}
	# Simulate traits
	trait <- sapply(seq_len(parameters$n.traits), function(x) phytools::fastBM(tree.lambda))
	colnames(trait) <- colnames(trait, do.NULL = FALSE, prefix = "t.")
	# Simulate correlation between traits
	if(parameters$dependent & parameters$n.traits>1){
		mc <- matrix(parameters$cor.traits, parameters$n.traits, parameters$n.traits)
		diag(mc) <- 1
		mchol <- chol(mc)
		trait <- trait%*%mchol
	}
	colnames(trait) <- colnames(trait, do.NULL = FALSE, prefix = "t.")
	# K statistic
	K.stat <- apply(trait, 2, function(x) picante::Kcalc(x = x, phy = tree, checkdata = FALSE))
	return(K.stat)
}

# Function to repeat the simulation with each parameter combination (K statistic)
# parameters = A data.frame with only one row and seven columns, the names of columns are all parameters setting to run the simulation
# n.rep = The number of repetitions
RUN.K <- function(parameters, n.rep){
	f <- function(x, parameters){
		RES <- get.K(parameters)
		return(RES)
	}
	RES <- t(sapply(as.list(1:n.rep), f, parameters = parameters, simplify = TRUE))
	return(RES)
}


# Function for data simulation, imputation by random value and to estimate NRMSE
# parameters = A data.frame with only one row and eight columns, the names of columns are all parameters setting to run the simulation
get.Random.Error <- function(parameters){
	# Simulate the phylogenetic tree
	tree <- geiger::drop.extinct(geiger::sim.bdtree(b = parameters$birth, d = parameters$death, stop = "taxa", n = parameters$n.spp, extinct = FALSE))
	tree$tip.label <- sprintf("Sp_%.4d", 1:length(tree$tip.label))	
	tree <- ape::makeNodeLabel(tree, method = "number", prefix = "Node")
	# Transform tree using Lambda or Grafen methods
	if(parameters$lambda<=1){
		tree.lambda <- geiger::rescale(tree, model = "lambda", parameters$lambda)
	} else{
		tree.lambda <- ape::compute.brlen(tree, method = "Grafen", power = parameters$lambda)
	}
	# Decomposing phylogenetic distance matrix derived from tree into a set of orthogonal vectors
	x <- PVR::PVRdecomp(tree)
	pvrs <- x@Eigen$vectors
	# Simulate traits
	trait <- sapply(seq_len(parameters$n.traits), function(x) phytools::fastBM(tree.lambda))
	colnames(trait) <- colnames(trait, do.NULL = FALSE, prefix = "t.")
	# Simulate correlation between traits
	if(parameters$dependent & parameters$n.traits>1){
		mc <- matrix(parameters$cor.traits, parameters$n.traits, parameters$n.traits)
		diag(mc) <- 1
		mchol <- chol(mc)
		trait <- trait%*%mchol
	}
	colnames(trait) <- colnames(trait, do.NULL = FALSE, prefix = "t.")
	# Introduce missing values
	traits.NA <- missForest::prodNA(trait, noNA = parameters$noNA)
	# Imputation random value
	n.NA.temp <- apply(is.na(traits.NA), 2, sum, na.rm = TRUE)
	range.temp <- apply(traits.NA, 2, range, na.rm = TRUE)
	traits.imp <- traits.NA
	for(i in which(n.NA.temp!=0)){
		traits.imp[is.na(traits.NA[,i]), i] <- runif(n = n.NA.temp[i], min = range.temp[1,i], max = range.temp[2,i])
	}
	# NRMSE
	nrmse.random <- missForest::nrmse(ximp = traits.imp, xmis = traits.NA, xtrue = trait)
	## Organize results
	results <- cbind(nrmse.random)
	colnames(results) <- c("Random")
	rownames(results) <- c("NRMSE")
	RES <- list()
	RES$Parameters <- parameters
	RES$Error <- results
	return(RES)
}

# Function to repeat the simulation with each parameter combination (random imputation)
# parameters = A data.frame with only one row and eight columns, the names of columns are all parameters setting to run the simulation
# n.rep = The number of repetition
# parallel = Number of parallel processes (default parallel = NULL)
RUN.random <- function(parameters, n.rep, parallel = NULL){
	f <- function(x,parameters){
		RES <- get.Random.Error(parameters)
		return(RES)
	}
	newClusters <- FALSE
	if(is.numeric(parallel)){
		parallel <- parallel::makeCluster(parallel, type = "PSOCK")
		newClusters <- TRUE
		parallel::clusterExport(parallel, "get.Random.Error")
	}
	if(!inherits(parallel, "cluster")){
		RES <- sapply(as.list(1:n.rep), f, parameters = parameters, simplify = FALSE)
	} else {
		RES <- parallel::parSapply(parallel, as.list(1:n.rep), f, parameters = parameters, simplify = FALSE)
	}
	if(newClusters) {
		parallel::stopCluster(parallel)
	}
	return(RES)
}
