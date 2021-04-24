## Using phylogenetic information to impute missing functional trait values in ecological databases

# Main functions

# Function to data simulation, imputation by missForest and to get the NRMSE
# parameters = A data.frame with only one row and eight columns, the names of columns are all parameters setting to run the simulation
get.missForest.Error <- function(parameters){
	# Simulate the phylogenetic tree
	tree <- geiger::drop.extinct(geiger::sim.bdtree(b = parameters$birth, d = parameters$death, stop = "taxa", n = parameters$n.spp, extinct = FALSE))
	tree$tip.label <- sprintf("Sp_%.4d", 1:length(tree$tip.label))	
	tree <- ape::makeNodeLabel(tree, method = "number", prefix = "Node")
	# Transform tree using Pagel lambda and Grafen rho methods
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
	# Imputation without PVR
	traits.imp <- missForest::missForest(traits.NA, maxiter = 15, ntree = 100, variablewise = FALSE)
	# NRMSE
	nrmse0 <- missForest::nrmse(ximp = traits.imp$ximp, xmis = traits.NA, xtrue = trait)
	## Imputation with PVRs
	nrmsePVR <- matrix(NA, 1, ncol(pvrs))
	for(i in 1:ncol(pvrs)){
		# Imputation
		traits.imp.pvr.temp <- missForest::missForest(cbind(traits.NA, pvrs[,1:i, drop = FALSE]), maxiter = 15, ntree = 100, variablewise = FALSE)
		# NRMSE
		nrmsePVR[1, i] <- missForest::nrmse(ximp = traits.imp.pvr.temp$ximp[, 1:ncol(trait)], xmis = traits.NA, xtrue = trait)
	}
	## Merge results
	results <- cbind(nrmse0, nrmsePVR)
	rownames(results) <- c("NRMSE")
	colnames(results) <- paste0("nPVR", 0:ncol(pvrs))
	RES <- list()
	RES$Parameters <- parameters
	RES$Error <- results
	return(RES)
}

# Function to repeat the simulation with each parameter combination (missForest imputation)
# parameters = A data.frame with only one row and eight columns, the names of columns are all parameters setting to run the simulation
# n.rep = The number of repetition
# parallel = Number of parallel processes (default parallel = NULL)
RUN.missForest <- function(parameters, n.rep, parallel = NULL){
	f <- function(x, parameters){
		RES <- get.missForest.Error(parameters)
		return(RES)
	}
	newClusters <- FALSE
	if(is.numeric(parallel)){
		parallel <- parallel::makeCluster(parallel, type = "PSOCK")
		newClusters <- TRUE
		parallel::clusterExport(parallel, "get.missForest.Error")
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

