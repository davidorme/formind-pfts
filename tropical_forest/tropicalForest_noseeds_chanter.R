source('../FORMIND_tree_models.R')

# Matching up the R functions to a simulation using the tropicalForests parameter file 
# with the following changes:
# - global seed ingrowth turned off
# - no seedlings permitted
# - background mortality turned off
# - stem diameter threshold for inclusion in the res file set to zero
# - max age set to 2000 to max out the climax PFT
# - growth curves switched to broadly similar shapes parameterised using 
#   the chanter model - easier to check parameterisation

# The pin file for this simulation contains one single initial seedling of each species, each
# in a separate patch

noseed_trees <- read.delim('tropicalForest_noseeds_chanter.res', skip=2)
noseed_trees <- split(noseed_trees, f=noseed_trees$Grp)
noseed_trees <- lapply(noseed_trees, function(x) x[1:which.max(x$H),])

names(noseed_trees) <- c('pioneer','intermed','climax')

# GEOMETRY check

PFTs <- list('pioneer' = list(Hmax = 20.0,
							  D_to_H = list(method='power', pars=list(h_0=23.87, h_1=0.49)),
							  D_to_Cd = list(method='linear', pars=list(cd_0=15.0)),
							  H_to_Clp = list(method='linear', pars=list(cl_0=0.25)),
							  Geom_to_B = list(method='geometric', pars=list(rho=0.4, sigma=0.7, 
												f_0=0.22, f_1=-0.22, form_factor=16)),
							  D_to_LAI = list(method='power', pars=list(l_0=2, l_1=0.0)),
							  D_to_Dinc = list(method='chanter', pars=list(alpha_0=0.42, alpha_1=2.1)),
							  GPP = list(pars=list(alpha=0.15, pmax=15.0, k=0.6, m=0.1, active=365)),
							  r_g = 0.2),
			 'intermed' = list(Hmax = 30.0,
								 D_to_H = list(method='power', pars=list(h_0=30.03, h_1=0.54)),
								 D_to_Cd = list(method='linear', pars=list(cd_0=15.0)),
								 H_to_Clp = list(method='linear', pars=list(cl_0=0.25)),
								 Geom_to_B = list(method='geometric', pars=list(rho=0.6, sigma=0.7, 
												  f_0=0.22, f_1=-0.22, form_factor=16)),
								 D_to_LAI = list(method='power', pars=list(l_0=2, l_1=0.0)),
								 D_to_Dinc = list(method='chanter', pars=list(alpha_0=0.074, alpha_1=1.5)),
								 GPP = list(pars=list(alpha=0.25, pmax=7.0, k=0.6, m=0.1, active=365)),
								 r_g = 0.2),
			 'climax'   = list(Hmax = 40.0,
								 D_to_H = list(method='power', pars=list(h_0=32.96, h_1=0.56)),
								 D_to_Cd = list(method='linear', pars=list(cd_0=15.0)),
								 H_to_Clp = list(method='linear', pars=list(cl_0=0.25)),
								 Geom_to_B = list(method='geometric', pars=list(rho=0.8, sigma=0.7, 
												  f_0=0.22, f_1=-0.22, form_factor=16)),
								 D_to_LAI = list(method='power', pars=list(l_0=2, l_1=0.0)),
								 D_to_Dinc = list(method='chanter', pars=list(alpha_0=0.0215, alpha_1=1.1)),
								 GPP = list(pars=list(alpha=0.35, pmax=4.0, k=0.6, m=0.1, active=365)),
								 r_g = 0.2))

# Find the appropriate maxima
PFTs <- lapply(PFTs, set_maxima)

# Overlay geometry onto the FORMIND res outputs
par(mfrow=c(5,2), mar=c(3,3,1,1), mgp=c(2,0.8,0))

plot_geometry_overlay(noseed_trees$pioneer, PFTs$pioneer)
plot_geometry_overlay(noseed_trees$intermed, PFTs$intermed)
plot_geometry_overlay(noseed_trees$climax, PFTs$climax)

# PRODUCTION CALCULATION

# DIAMETER INCREMENT 
# Needs to allow for the offset - the res file puts the biomass increment
# from Year N in the row for Year N+1 (how much has this grown since we last 
# saw it), whereas the function predicts the growth for this year.


# Now, does the D in N+1 simply step up by DInc? 
# Note that in the res files, DInc is the gain _from the 
# previous year_ so D_N+1 = D_N + DInc_N?

par(mfrow=c(3,2))
for(df in noseed_trees){

	# Get the offsets right. Final Year has an unknown gain from the next year
	N <- nrow(df)
	df$D2 <- c(df$D[1], (df$D[-N] + df$DInc[-1]))
	
	# Look at proportion difference
	plot(D2 ~ D, data=df)
	abline(0,1, col='red')
	plot(density(df$D - df$D2))
}

# Replicate the calculation of growth using the R code

for(pft in names(PFTs)){
	tmp <-  D_to_Dinc(noseed_trees[[pft]]$D, PFTs[[pft]])
	# offset DIncR to the following year
	noseed_trees[[pft]]$DIncR <- c(0, tmp[-length(tmp)])
}

par(mfrow=c(3,2))
for(df in noseed_trees){
	plot(DInc ~ D, data=df, type='l', lwd=3, col='grey')
	lines(DIncR ~ D, data=df, col='red')
	with(df, hist(DInc - DIncR))
}

# The pioneer and climax are spot on, to within rounding error given
# reported decimal places, but oddly intermed is slightly out. 

# BIOMASS INCREMENT

# Mechanism check: BT_N+1 = BT_N + BInc_N * (1 - growth respiration)

par(mfrow=c(3,2))
for(df in noseed_trees){
	
	# - Omit the last year
	N <- nrow(df)
	df$BT2 <- c(df$BT[1], df$BT[-N] + df$BInc[-1])
	
	# Look at correspondance and proportion difference 
	# - precise to limits of rounding
	plot(BT ~ BT2, data=df)
	abline(0,1, col='red')
	plot(density(df$BT - df$BT2))
}

# Replicate this using the R code.
# Get new height, given D increment, get new Biomass: Geom_to_B(D + DInc, H + HInc) - Geom_to_B(D, H)
for(pft in names(PFTs)){
	
	N <- nrow(noseed_trees[[pft]])
	
	# calculate the Biomass from the geometry
	BT <- Geom_to_B(noseed_trees[[pft]]$D, noseed_trees[[pft]]$H, PFTs[[pft]])
	#print(cor(BT, noseed_trees[[pft]]$BT), digits=20)
	
	# Get the new geometry
	DNew <- noseed_trees[[pft]]$D[-N] + noseed_trees[[pft]]$DIncR[-1]
	HNew <- D_to_H(DNew, PFTs[[pft]])
	BTNew <- Geom_to_B(DNew, HNew, PFTs[[pft]])
	# print(cor(BTNew, noseed_trees[[pft]]$BT[-1]), digits=20)
	
	noseed_trees[[pft]]$BIncR <- c(0,  BTNew - BT[-N])
		
	# The respiration growth is then the GPP that would be needed
	# to attain the BT Increment. This is parameterised as a proportion
	# of total growth spend, so needs to be converted to a proportion of
	# the actual growth.
	r_g = 0.2
	noseed_trees[[pft]]$R_GrowthR <- noseed_trees[[pft]]$BIncR * (r_g / (1 - r_g))
	
}

# Plot those against the FORMIND BInc - intermed is off
par(mfrow=c(3,2))
for(df in noseed_trees){
	plot(BInc ~ D, data=df, type='l', lwd=3, col='grey')
	lines(BIncR ~ D, data=df, col='red')
	with(df, hist(BInc - BIncR))
}

par(mfrow=c(3,2))
for(df in noseed_trees){
	plot(R_Growth ~ D, data=df, type='l', lwd=3, col='grey')
	lines(R_GrowthR ~ D, data=df, col='red')
	with(df, hist(R_Growth - R_GrowthR))
}

# Compare the growth through time from each mechanism

pdf('Tropical_forest_growth.pdf', paper='a4', height=8, width=8)
	par(mfrow=c(3,2), mar=c(3,3,1,1), mgp=c(2,0.8,0))
	
	for(pft in names(PFTs)){
		
		res <- noseed_trees[[pft]]
		plot_production(res, main= sprintf('Formind: %s', pft))
		
		sim <- growth_simulation(PFTs[[pft]], Dinit=0.001, n_years=max(res$AGE+1))
		plot_production(sim, main= sprintf('R: %s', pft))
	}
dev.off()

