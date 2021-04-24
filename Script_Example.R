## Using phylogenetic information to impute missing functional trait values in ecological databases
# Simplified example to imputation procedure in real incomplete dataset
# See cited references in the main text, e.g. Penone et al. 2014

## Load packages
require(geiger)
require(PVR)
require(phytools)
require(missForest)

## Data simulation for this example
# Phylogenetic tree
tree <- geiger::drop.extinct(geiger::sim.bdtree(b = 0.1, 
												d = 0, 
												stop = "taxa", 
												n = 50, 
												extinct = FALSE))
# Trait dataset
trait <- sapply(seq_len(2), function(x) phytools::fastBM(tree))
colnames(trait) <- colnames(trait, do.NULL = FALSE, prefix = "t.")
# Introduce missing values
traits.NA <- missForest::prodNA(trait, noNA = 0.25)

## Imputation
# Decomposing phylogenetic distance matrix derived from tree into a set of orthogonal vectors
x <- PVR::PVRdecomp(tree)
# Extract the PVRs
pvrs <- x@Eigen$vectors

# Combine traits and PVRs
traits.pvrs <- cbind(traits.NA, pvrs)

# Imputation using missForest (note that this function have other arguments, see details in Stekhoven and Buehlmann 2012)
RF.imp <- missForest::missForest(traits.pvrs, maxiter = 15, ntree = 100, variablewise = FALSE)

# OOBerror
RF.imp$OOBerror

# Imputed dataset
traits.imp <- RF.imp$ximp[, seq_len(ncol(traits.NA)), drop = FALSE]
traits.imp
