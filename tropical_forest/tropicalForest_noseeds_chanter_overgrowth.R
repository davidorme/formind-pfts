source('../FORMIND_tree_models.R')

# Plotting the res file outputs for the overgrown chanter parameterisation

# with the following changes:
# - global seed ingrowth turned off
# - no seedlings permitted
# - background mortality turned off
# - stem diameter threshold for inclusion in the res file set to zero
# - max age set to 2000 to max out the climax PFT
# - growth curves switched to _wildly overclocked_ models parameterised using 
#   the chanter model - easier to check parameterisation

# The pin file for this simulation contains one single initial seedling of each species, each
# in a separate patch

noseed_trees <- read.delim('tropicalForest_noseeds_chanter_overgrowth.res', skip=2)
noseed_trees <- split(noseed_trees, f=noseed_trees$Grp)
noseed_trees <- lapply(noseed_trees, function(x) x[1:which.max(x$H),])

names(noseed_trees) <- c('pioneer','intermed','climax')


# PRODUCTION plots
chanter_params = list('pioneer'= list(dmax=max(noseed_trees[['pioneer']]$D), alpha_0=1, alpha_1=2.1),
					  'intermed'= list(dmax=max(noseed_trees[['intermed']]$D), alpha_0=0.4, alpha_1=1.5),
					  'climax'= list(dmax=max(noseed_trees[['climax']]$D), alpha_0=0.1, alpha_1=1.1))

pdf('Tropical_forest_growth_overgrowth.pdf', paper='a4', height=11, width=8)
	

	par(mfrow=c(3,2), mar=c(3,3,1,1), mgp=c(2,0.8,0))
	for(pft in names(noseed_trees)){
		
		# data
		res <- noseed_trees[[pft]]
		pars <- chanter_params[[pft]]
		
		# plot the DInc curves
		D <- seq(0, pars$dmax, length=101)
		DInc <- with(pars, alpha_0 * D * ( 1 - D/dmax) * exp(-alpha_1 * D))
		plot(DInc ~ D, type='l')
		lines(DInc[-1] ~ D[-nrow(res)], data=res, col='red')
		
		# plot the production curves
		plot_production(res, main= sprintf('Formind production: %s', pft))
		
	}
dev.off()

