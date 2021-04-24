## Using phylogenetic information to impute missing functional trait values in ecological databases
# Code for main results, tables S4-S9 and figure 2

## Load packages
require(ape)
require(geiger)
require(PVR)
require(phytools)
require(missForest)
require(parallel)

## Load functions
source("Main_functions.R")

## Set all parameters
# In this example trees are set with 100 species and dataset with 3, 5 and 10 traits (tables S7-S9)
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

## Run the simulation and save the NRMSE in .csv file

# Create directory to save the results
dir.create("Results")

# Set main parameters to the simulation
# This can take a long time to finish

# Set main parameters to the simulation
n.rep <- 1000 # Number of repetition
parallel <- NULL # Number of parallel processes, it depends on the computer that is running the analysis

for(i in 1:nrow(PAR)){
	print(i)
	parameters <- PAR[i,]
	# Run missForest imputation to get the NRMSE
	results.temp <- RUN.missForest(parameters, n.rep = n.rep, parallel = parallel)
	# Organize the results
	results <- sapply(results.temp, function(x) x$Error[1, 1:(parameters$n.spp-1), drop = FALSE], simplify = FALSE)
	results <- do.call(rbind, results)
	rownames(results) <- paste0("list.", seq_len(n.rep), ".nrmse")
	# Base file name
	base <- paste(names(parameters), parameters, sep = "_", collapse = "_")
	# Write the .csv file
	write.csv(results, paste0("./Results/RES_", sprintf("%.4d", i), "_nrmse_", base, ".csv"))
}


## Organize the results

# Small function to organize the results
f1 <- function(x){
	res<- c(x[1], # NRMSE without PVR
			min(x[-1]), # Minimum NRMSE with PVR
			which(x[-1]==min(x[-1])), # Number of PVR in the minimum NRMSE
			ifelse(x[1]<min(x[-1]), 1, 0)) # 1 if NRMSE without PVR is smaller that minimum NRMSE with PVR
	names(res) <- c("nrmse.without.pvr", 
					"nrmse.min.with.pvr",
					"n.pvr.min",
					"smaller.without.pvr")
	return(res)
}

# Small function to get summary statistics
f2 <- function(x){
	r.mean.temp <- apply(x, 2, mean)
	r.sd.temp <- apply(x, 2, sd)
	r.median.temp <- apply(x, 2, median)
	r.sum.temp <- apply(x, 2, sum)
	res <- c(r.mean.temp[1], r.sd.temp[1], # Mean and standard deviation to NRMSE without PVR
			 r.mean.temp[2], r.sd.temp[2], # Mean and standard deviation to NRMSE with PVR
			 r.median.temp[3], # Median to number of PVR in the minimum NRMSE
			 r.sum.temp[4]/nrow(x)) # Proportion simulation with smaller NRMSE without PVR
	names(res) <- c("mean.nrmse", 
					"sd.nrmse", 
					"mean.nrmse.pvr", 
					"sd.nrmse.pvr",
					"median.n.pvr",
					"prop.sim.smaller.without.pvr")
	return(res)
}


# Organize the results of NRMSE
RESULTS.nrmse <- matrix(NA, nrow(PAR), ncol(PAR)+6)
for(i in 1:nrow(PAR)){
	print(i)
	parameters <- PAR[i,]
	# Base file name
	base <- paste(names(parameters), parameters, sep = "_", collapse = "_")
	# Read .csv with the results
	res.nrmse <- read.csv(paste0("./Results/RES_", sprintf("%.4d", i), "_nrmse_", base, ".csv"), row.names = 1)
	# Organize the results
	res.temp.nrmse <- t(apply(res.nrmse, 1, f1))
	# Summary statistic of results
	res.temp.nrmse.summary <- f2(res.temp.nrmse)
	RESULTS.nrmse[i, ] <- unlist(c(parameters, res.temp.nrmse.summary))
	# Set names to RESULTS.nrmse object
	if(i==1){
		colnames(RESULTS.nrmse) <- names(unlist(c(parameters, res.temp.nrmse.summary)))
	}
}

## Final results to NRMSE
RESULTS.nrmse <- as.data.frame(RESULTS.nrmse)
RESULTS.nrmse


## Plot the results

# Load packages to plot
require(magrittr)
require(dplyr)
require(plotrix)


## Plot function
# This function include a lot of options and arguments to produce the figure 2
f.plot.2 <- function(results.temp, legend = FALSE, xlab = FALSE, ylab = FALSE, lim.y = 1.25, opt.adj = 0.05, leg.adj = 1.3){
	x <- as.numeric(as.factor(results.temp$noNA))
	adj <- as.numeric(as.factor(results.temp$lambda))*0.1
	plot(x, type = "n", xaxt = "n", 
		 ylab = "", xlab = "",
		 xlim = c(0.5, max(x)+0.5),
		 ylim = c(0,lim.y), cex.axis = 1.2)
	grid(nx = NA, ny = NULL, col = "lightgray", lty = 1)
	if(xlab){
		title(xlab = "Proportion of missing data", cex.lab = 1.5)  
	}
	if(ylab){
		title(ylab = "NRMSE", cex.lab = 1.3)
	}
	axis(1, seq_len(length(levels(as.factor(results.temp$noNA)))), levels(as.factor(results.temp$noNA)), cex.axis = 1.2)
	palette(c("#ffa06d", "#ff6e40", "#dd2c00", "#a30000"))
	plotCI(x-0.4+adj, results.temp$mean.nrmse, uiw = results.temp$sd.nrmse, 
		   col = as.factor(results.temp$lambda), pch = 19, add = TRUE)
	palette(c("#63ccff", "#039be5", "#01579b", "#002f6c"))
	plotCI(x+0+adj, results.temp$mean.nrmse.pvr, uiw = results.temp$sd.nrmse.pvr, 
		   col = as.factor(results.temp$lambda), pch = 15, add = TRUE)
	if(legend){
		palette(c("#63ccff", "#039be5", "#01579b", "#002f6c"))
		legend(4.8, leg.adj, title = "Phylogenetic \n signal", 
			   legend = c(expression(lambda == 0.2), expression(lambda == 0.6), expression(lambda == "1.0"), expression(rho == 1.2)), 
			   pch = c(15), col = palette(), bty = "o", inset = 0.05, pt.cex = 1.3, cex = 1.3,
			   bg = "white", box.col = "white")
		palette(c("#ffa06d", "#ff6e40", "#dd2c00", "#a30000"))
		legend(4.7, leg.adj, title = "", legend = c("    ","","",""), 
			   pch = c(19), col = palette(), bty = "n", inset = 0.05, pt.cex = 1.3, cex = 1.3)
		legend(2, leg.adj, title = "missForest \n method", legend = c("without phylogeny", "with phylogeny"), 
			   pch = c(19, 15), col = "grey40", bty = "o", inset = 0.05, pt.cex = 1.3, cex = 1.3, title.adj = 0.3,
			   bg = "white", box.col = "white")
	}
	xt <- x+0+adj
	yt <- results.temp$mean.nrmse.pvr+results.temp$sd.nrmse.pvr+opt.adj
	sw   <- strwidth(round(results.temp$median.n.pvr))
	sh   <- strheight(round(results.temp$median.n.pvr))
	frsz <- 0.01
	rect(xt - sw/2 - frsz,
		 yt - sh/2 - frsz,
		 xt + sw/2 + frsz,
		 yt + sh/2 + frsz,
		 col = "grey90", border = NA)
	text(x+0+adj, results.temp$mean.nrmse.pvr+results.temp$sd.nrmse.pvr, 
		 labels = round(results.temp$median.n.pvr),
		 pos = 3, cex = 1)
	if(unique(results.temp$dependent)==0){
		mtext("Correlation = Independent", line = 0, adj = 0, font = 2, cex = 0.8)
	}
	if(unique(results.temp$dependent)==1){
		mtext(paste0("Correlation = ", unique(results.temp$cor.traits)), line = 0, adj = 0, font = 2, cex = 0.8)
	}
}


# Filter the results to plot

# Set number of species and number of traits to filter the results
nTraits <- 10
nSpp <- 100

names(RESULTS.nrmse)

results.temp1 <- RESULTS.nrmse %>% filter(n.traits == nTraits,
										  n.spp == nSpp,
										  dependent == 1,
										  cor.traits == 0.10)
results.temp2 <- RESULTS.nrmse %>% filter(n.traits == nTraits,
										  n.spp == nSpp,
										  dependent == 1,
										  cor.traits == 0.25)
results.temp3 <- RESULTS.nrmse %>% filter(n.traits == nTraits,
										  n.spp == nSpp,
										  dependent == 1,
										  cor.traits == 0.50)
results.temp4 <- RESULTS.nrmse %>% filter(n.traits == nTraits,
										  n.spp == nSpp,
										  dependent == 1,
										  cor.traits == 0.75)
results.temp5 <- RESULTS.nrmse %>% filter(n.traits == nTraits,
										  n.spp == nSpp,
										  dependent == 1,
										  cor.traits == 0.90)
results.temp6 <- RESULTS.nrmse %>% filter(n.traits == nTraits,
										  n.spp == nSpp,
										  dependent == 0,
										  cor.traits == 0)

# Export the plot in .pdf
pdf(paste0("Res_nrmse_nSpp_", nSpp, "_nTraits_", sprintf("%.2d", nTraits), ".pdf"), height = 12, width = 8)
par(mfrow = c(3,2), mar = c(4.5, 4, 0.5, 0.5), oma = c(1.1, 1, 2.5, 1), las = 1)
f.plot.2(results.temp6, ylab = TRUE, lim.y = 1.25, opt.adj = 0.045)
f.plot.2(results.temp1, lim.y = 1.25, opt.adj = 0.045)
f.plot.2(results.temp2, ylab = TRUE, lim.y = 1.25, opt.adj = 0.045)
f.plot.2(results.temp3, lim.y = 1.25, opt.adj = 0.045)
f.plot.2(results.temp4, ylab = TRUE, xlab = TRUE, lim.y = 1.25, opt.adj = 0.045)
f.plot.2(results.temp5, legend = TRUE, xlab = TRUE, lim.y = 1.25, opt.adj = 0.045, leg.adj = 1.21)
dev.off()

