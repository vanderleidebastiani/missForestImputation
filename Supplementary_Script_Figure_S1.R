## Using phylogenetic information to impute missing functional trait values in ecological databases
# Code for results of the figure S1

## Load packages
require(ape)
require(geiger)
require(PVR)

## Load functions
source("Supplementary_functions.R")

## Set all parameters
# Figure S1
PAR <- expand.grid(birth = c(0.1), # Speciation rate in the phylogenetic tree
				   death = c(0), # Extinction rate in the phylogenetic tree
				   n.spp = c(100))  # Number of species in the phylogenetic tree
PAR

## Run the simulation 

# Set main parameters for the simulation
n.rep <- 1000 # Number of repetitions

RESULTS.pvr <- RUN.PVR(PAR, n.rep = n.rep)


## Organize the results
RESULTS.pvr <- do.call(rbind, RESULTS.pvr)

# Summary statistic of results
results.mean <- apply(results, 2, mean)
names(results.mean) <- paste0("mean.", names(results.mean))
results.sd <- apply(results, 2, sd)
names(results.sd) <- paste0("sd.", names(results.sd))

## Final results
results.mean
results.sd

## Plot

# Load packages to plot
require(plotrix)

# Select only the first 50 PVRs
pvr.mean <- results.mean[1:50]
pvr.sd <- results.sd[1:50] 

# Export the plot in .pdf
pdf("Res_PVRs_100spp.pdf")
par(mar = c(5.2, 4.3, 2, 2))
plot(pvr.mean, type = "n", xaxt = "n", 
	 xlab = "", ylab = "",
	 xlim = c(1, 50),
	 ylim = c(0, 1),
	 cex.axis = 1.2)
plotCI(seq_len(length(pvr.mean)), pvr.mean, uiw = pvr.sd, 
	   col = "#a30000", pch = 19, add = T)
title(xlab = "Eigenvector number", cex.lab = 1.5)
title(ylab = "Cumulative Eigenvalue", cex.lab = 1.3)
axis(1, c(1,10,20,30,40,50), cex.axis = 1.2)
dev.off()

