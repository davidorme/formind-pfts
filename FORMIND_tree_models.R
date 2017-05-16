## FUNCTIONS TO EXPLORE THE FORMIND MODELS FOR TREE GEOMETRY AND PRODUCTION.

## GEOMETRY SCALING FUNCTIONS. These are the geometric relationships that an individual will
## follow but their annual progress along the curves is set by the productivity 

set_maxima <- function(settings){
	
	# Function to update a settings object to set the maxima
	# - Dmax can be set directly or via the maximum attainable height
	# - Hmax can be set directly or via the maximum attainable diameter
	
	# what is provided
	hd <- c(is.null(settings$Hmax), is.null(settings$Dmax))
	hd <- paste(ifelse(hd, c('-','-'), c('h', 'd')), collapse='')
	
	if(hd == '--'){
		stop('Settings must contain one of Hmax or Dmax')
	} else if(hd == 'hd'){
		warning('Both Hmax and Dmax already set.')
	} else if(hd == '-d'){
		# set temporary Hmax as it is required in D_to_H
		settings$Hmax <- Inf
		settings$Hmax <- D_to_H(settings$Dmax, settings)
	} else if(hd == 'h-'){
		# Find D at which Hmax - D_to_H = 0 
		fn <- function(d, local_set){
			# store the target value and then _locally_ alter Dmax and make Hmax
			# infinite to provide a monotonic function for uniroot
			target <- local_set$Hmax
			local_set$Dmax <- Inf
			local_set$Hmax <- Inf
			return(target - D_to_H(d, local_set))
		}
		settings$Dmax = uniroot(fn, local_set=settings, interval=c(0,1e5))$root
	}
	
	return(settings)
}

D_to_H <- function(D, settings, verbose=TRUE){
	
	# Stem diameter (D) to height (H) functions
	
	# Warn about D > Dmax (unless we're running in silent mode
	# when uniroot might be being used.)
	if(any(D > settings$Dmax) & verbose) warning('Some D values greater than Dmax')
	
	H <- with(settings$D_to_H,
			with(pars, switch(method,
				'power' =  h_0 * D ^ h_1,
				'saturation' = D/((1/h_0)+(D/h_1)),
				'polynomial' = h_0 + h_1 * D + h_2 * D^2)))
	
	return(pmin(H, settings$Hmax))
}

D_to_Cd <- function(D, settings){
	
	# Stem diameter (D) to Crown diameter (CD)
	
	Cd <- with(settings$D_to_Cd,
			with(pars, switch(method,
				'exponential_1' = D * (cd_0 + cd_1 * exp(-cd_2 * D)),
				'exponential_2' = cd_0 * D + cd_1 * exp(-cd_2 * D),
				'polynomial' = cd_0 + cd_1 * D + cd_2 * D^2 + cd_3 * D^3,
				'linear' = cd_0 * D,
				'saturation' = D/((1/cd_0)+(D/cd_1)),
				'power' = cd_0 * D ^ cd_1 -cd_2)))
	
	return(Cd)
}

H_to_Clp <- function(H, settings){
	
	# Height (H) to Crown length proportion (CLP) functions
	
	Clp <- with(settings$H_to_Clp,
			with(pars, switch(method,
				'linear' = cl_0,
				'saturation' = (- ((cl_0 * H * cl_1)/(cl_0 * H + cl_1))),
				'polynomial' = (cl_0 + cl_1 * H + cl_2 * H^2))))
	
	return(Clp)
}

Cd_to_Ca <- function(Cd){
	
	# Crown diameter (Cd) to Crown area (Ca)
	
	return(pi/4 * Cd^2)
}

D_to_LAI <- function(D, settings){
	
	# Stem diameter (D) to LAI functions
	
	LAI <- with(settings$D_to_LAI, 
			with(pars, switch(method,
				'linear' = l_0 + l_1 * (D/100),
				'power' = l_0 * D ^ l_1)))
	
	return(LAI)
}

Geom_to_B <- function(D, H, settings){
	
	# Geometry to Biomass (B) functions. more moving parts in this one
	B <- with(settings$Geom_to_B,
			with(pars, switch(method,
				'power' = b_0 * D ^ b_1,
				'logarithmic' = {
					Dt <- log(D) - b_2
					B <- exp(b_0 * (Dt) * ((2 * b_1 + Dt)/(b_1 + Dt)))
				},
				'geometric' = {
					if(form_factor == 15){
						f <- f_0 * exp(f_1 * D ^ f_2)
					} else if(form_factor == 16) {
						f <- f_0 * D ^ f_1
					} else {
						stop('Unknown form factor equation')
					}
			
					B <- (pi/4) * D^2 * H * f * (rho/sigma)
				})))
	
	return(B)
}

geometry_curves <- function(settings, N=101){
	
	# Ensemble function to create a dataframe of the geometry
	# over the range of stem diameters
	
	D <- seq(0, settings$Dmax, len=N)
	
	H <- D_to_H(D, settings)
	CD <- D_to_Cd(D, settings)
	CLP <- H_to_Clp(H, settings)
	CA <- Cd_to_Ca(CD)
	BT <- Geom_to_B(D, H, settings)
	LAI <- D_to_LAI(D, settings)
	
	return(data.frame(D=D, H=H, CD=CD, CLP=CLP, CA=CA, BT=BT, LAI=LAI))
}

## PRODUCTIVITY FUNCTIONS. Gross primary productivity and diameter increment functions.

GPP <- function(LAI, CA, settings, irradiance=700, day_length=12, phi_T=1.0, phi_W=1.0){
	
	# Function to calculate the annual gross primary productivity 
	# from the incident light on a tree given the LAI and 
	# crown area and some PFT specific parameters:
	# - alpha: initial slope of light response curve
	# - pmax: maximum leaf gross photosynthetic rate
	# - k: light extinction coefficient
	# - m: transmission coefficient
	# - active: number of days per year actively growing
	
	# Interim photosynthesis is modified to GPP via two proportions
	# that represent water (phi_W) or temperature (phi_T) limitation
	
	# gross photosynthetic per second per m2 of canopy
	p_ind = with(settings$GPP$pars, (pmax/k) * log((alpha * k * irradiance + pmax*(1-m)) / 
									(alpha * k * irradiance * exp(-k*LAI)+ pmax*(1-m))))
	
	# - convert to whole canopy area per year
	p_ind = p_ind * CA * 3600 * day_length * settings$GPP$pars$active
	# - moles of photons to moles of CO2, grams of CO2, tonnes of ODM 
	p_ind = p_ind * 0.63 * 44 * 10^-12
	
	# Enforce any light or water limitations
	p_ind = p_ind * phi_W * phi_T
	
	return(p_ind)
}

D_to_Dinc <- function(D, settings){
	
	# Function to calculate the (maximal) diameter increment for a stem
	# given its current diameter D and growth parameterisation
	
	Dinc <- with(settings$D_to_Dinc,
				with(pars, switch(method,
					'polynomial' = {
						# calculate the internal variables needed
						dmax <- settings$Dmax
						dgmax <- (dpropgmax * dmax)
						x0 <- gmax - gstart * gmax
						x1 <- 2 * dgmax * (dmin - dmax) - dmin^2 + dmax^2
						x2 <- 2 * dgmax * (dmin - dgmax) - dmin^2 + dgmax^2
						x3 <- dgmax^4 * (dmax - dmin)
						x4 <- 2 * dgmax^3 * (dmin^2 - dmax^2)
						x5 <- dgmax^2 * (5*dmin^3 + 3 * dmin * dmax^2 - 3 * dmax * dmin^2 + dmax^3)
						x6 <- 2 * dgmax * (dmax * dmin^3 - dmin * dmax^3)
						x7 <- dmax^3 * dmin^2 - dmax^2 *dmin^3 + dmin^4 - dmin^5
						x8 <- 3 * dmin * dgmax^2 - 2 * dgmax^3 -dmin^3
			
						a3 <- (x0 * x1 + (gstart - gend) * gmax *x2) / (x3 + x4 + x5 + x6 + x7)
						a2 <- (x0 - a3 * x8) / (2 * dmin * dgmax - dgmax^2 - dmin^2)
						a1 <- -3 * a3 * dgmax^2 - 2 * a2 * dgmax
						a0 <- gstart * gmax - a3 * dmin^3 - a2 * dmin^2 - a1 * dmin
			
						a0 + a1 * D + a2 * D^2 + a3 * D^3},
					'weibull' = {alpha_0 * alpha_1 * alpha_2 * (alpha_1 * D) ^ (alpha_2 - 1) * exp(-(alpha_1 * D) ^ alpha_2)},
					'chanter' = {alpha_0 * D * ( 1 - D/settings$Dmax) * exp(-alpha_1 * D)})))
	
	return(Dinc)
}

## SIMULATION FUNCTION. Simple simulation under constant illumination and day length,
## although could easily be extended to take vectors of irradiance etc.

growth_simulation <- function(settings, Dinit=0.01, n_years=200, truncate_overgrowth=FALSE,
							  irradiance=700, day_length=12, phi_T=1.0, phi_W=1.0){
	
	# The growth and production outputs could all be calculated after the main loop
	# except when the growth function overshoots the available biomass, when
	# DInc needs to be recalculated to fit within bounds. It isn't sensible
	# to use a parameterisation that requires this, but for verifying against
	# FORMIND outputs, the behaviour is duplicated.
	
	# initialise the variable vectors
	Age <- 0:n_years
	D <- H <- CD <- CLP <- CA <- LAI <- BT <- PB <- DInc <- numeric(n_years+1)
	D[1] <- Dinit
	
	# Growth and flag for overgrowth
	overgrowth <- FALSE
	BInc <- R_Growth <- R_Main <- R <- numeric(n_years+1)
	
	# run the simulation.
	for(yr in seq_along(Age)){
		
		# fill in the rest of the geometry for this year
		H[yr] <- D_to_H(D[yr], settings, verbose=FALSE)
		CD[yr] <- D_to_Cd(D[yr], settings)
		CLP[yr] <- H_to_Clp(H[yr], settings)
		CA[yr] <- Cd_to_Ca(CD[yr])
		BT[yr] <- Geom_to_B(D[yr], H[yr], settings)
		LAI[yr] <- D_to_LAI(D[yr], settings)
		
		# now growth & production
		PB[yr] <- GPP(LAI[yr], CA[yr], settings, irradiance=irradiance, 
					  day_length=day_length, phi_T=phi_T, phi_W=phi_W)
		
		# Work out potential growth
		potDInc <- D_to_Dinc(D[yr], settings)
		potH <- D_to_H(D[yr] + potDInc, settings, verbose=FALSE)
		potBInc <- Geom_to_B(D[yr] + potDInc, potH, settings) - BT[yr]
		potR_Growth <- potBInc * (settings$r_g / (1 - settings$r_g))
		
		if((potR_Growth + potBInc) > PB[yr]){
			overgrowth <- TRUE
			if(truncate_overgrowth){
				# So we need to find a diameter increment that yields
				# a biomass increment that soaks up the available PB
				fn <- function(dInc, currentD, currentBT, targetPB, settings){
					
					d <- dInc + currentD
					h <- D_to_H(d, settings, verbose=FALSE)
					bt <- Geom_to_B(d, h, settings)
					binc <- bt - currentBT
					pb <- binc / (1 - settings$r_g)
					
					return(targetPB - pb)
				}
				
				# find the value where R_Growth and BInc soak up all of PB
				potDInc <- uniroot(fn, interval=c(0, potDInc), currentD=D[yr], currentBT=BT[yr],
								   targetPB=PB[yr], settings=settings)$root
				
				# recalculate the rest
				potH <- D_to_H(D[yr] + potDInc, settings)
				potBInc <- Geom_to_B(D[yr] + potDInc, potH, settings) - BT[yr]
				potR_Growth <- potBInc * (settings$r_g / (1 - settings$r_g))
				
			}
		}
		
		# start diameter for next year
		if(yr < (n_years+1)){
			D[yr + 1] <- D[yr] + potDInc
			DInc[yr + 1] <- potDInc
			BInc[yr + 1] <- potBInc
			R_Growth[yr + 1] <- potR_Growth
			R_Main[yr + 1] <- PB[yr] - (potBInc + potR_Growth)
			R[yr + 1] <- R_Growth[yr + 1] + R_Main[yr + 1]
		}
	}
	
	if(overgrowth) warning('Growth exceeded production!')

	return(data.frame(AGE=Age, D=D, H=H, CD=CD, CLP=CLP, CA=CA, BT=BT, LAI=LAI, 
					  PB=PB, DInc=DInc, BInc=BInc, R=R, R_Growth=R_Growth, R_Main=R_Main))
}

# # Older version that doesn't have the ability to constrain overgrowth
# # to actual production but is going to be faster when growth isn't 
# # constrained by production
# 
# growth_simulation <- function(settings, Dinit=0.01, n_years=200,
# 							  irradiance=700, day_length=12, phi_T=1.0, phi_W=1.0){
#
# 	# initialise the variable vectors
# 	Age <- seq_len(n_years)
# 	D <- H <- CD <- CLP <- CA <- LAI <- BT <- PB <- DInc <- numeric(n_years)
# 	D[1] <- Dinit
#
# 	# run the simulation.
# 	for(yr in Age){
#
# 		# fill in the rest of the geometry for this year
# 		H[yr] <- D_to_H(D[yr], settings)
# 		CD[yr] <- D_to_Cd(D[yr], settings)
# 		CLP[yr] <- H_to_Clp(H[yr], settings)
# 		CA[yr] <- Cd_to_Ca(CD[yr])
# 		BT[yr] <- Geom_to_B(D[yr], H[yr], settings)
# 		LAI[yr] <- D_to_LAI(D[yr], settings)
#
# 		# production
# 		PB[yr] <- GPP(LAI[yr], CA[yr], settings, irradiance=irradiance,
# 					  day_length=day_length, phi_T=phi_T, phi_W=phi_W)
#
# 		# start diameter for next year
# 		if(yr < n_years){
# 			D[yr + 1] <- D[yr] + D_to_Dinc(D[yr], settings)
# 		}
# 	}
#
# 	# Fill in the production
# 	DInc <- c(0, diff(D))
# 	BInc <- c(0, diff(BT))
# 	R_Growth <- BInc * (settings$r_g / (1 - settings$r_g))
#
# 	# Handle problems with growth exceeding production
# 	if(any((R_Growth + BInc) > PB)){
# 		warning('Growth parameterisation problem!')
# 	}
#
# 	# Calculate the remainders
# 	R_Main <- c(0, (PB[-n_years] - BInc[-1] - R_Growth[-1]))
# 	R <- R_Main + R_Growth
#
# 	return(data.frame(AGE=Age, D=D, H=H, CD=CD, CLP=CLP, CA=CA, BT=BT, LAI=LAI,
# 					  PB=PB, DInc=DInc, BInc=BInc, R=R, R_Growth=R_Growth, R_Main=R_Main))
# }


## PLOTTING FUNCTIONS. All will take either a chunk of a RES file for a single tree
## or the output of the growth_simulation or geometry_curves.

plot_production <- function(data, ...){
	
	# Function taking the output of the res file for a single tree
	# or the output of growth_dataulation and plotting the production
	# partitioning through time.
	
	plot(PB ~ AGE, data=data, type='l', lwd=3, col='grey', ...)
	lines(BInc ~ AGE, data=data, col='darkgreen', lty=2)
	lines(R ~ AGE, data=data, col='black', lty=3)
	lines(R_Main ~ AGE, data=data, col='red', lty=4)
	lines(R_Growth ~ AGE, data=data, col='blue', lty=5)

	legend('topleft', bty='n', lwd=c(3,1,1,1,1), lty=1:5,
			col=c('grey','darkgreen','black','red','blue'), 
			legend=c('Total biomass production (PB = BInc + R)', 'Growth (BInc)', 
					'Total respiration (R = Rm + Rg)', 'Maintenance respiration (Rm)', 
					'Growth respiration (Rg)'))
}

plot_geometry <- function(data, focal_row=0){
	
	# Function to plot out the six geometry variables from
	# either a geometry curves data frame or from a growth 
	# simulation or from a FORMIND simulation loaded from a 
	# res file.
	
	# The res files don't include canopy area
	if(is.null(data$CA)) data$CA <- Cd_to_Ca(data$CD)
	
	with(data,{
		
		Hmax <- max(H)
		Dmax <- max(D)
		
		plot(H ~ D, xlab='Stem diameter [m]', ylab='Height [m]', type='l')
		abline(h=Hmax, v=Dmax, col='grey', lty=2)
		if(focal_row) points(D[focal_row], H[focal_row], col='red')
	
		plot(CLP ~ H, xlab='Stem height [m]', ylab='Crown length proportion [-]', type='l')
		abline(v=Hmax, col='grey', lty=2)
		if(focal_row) points(H[focal_row], CLP[focal_row], col='red')
	
		plot(CD ~ D, xlab='Stem diameter [m]', ylab='Crown diameter [m]', type='l')
		abline(v=Dmax, col='grey', lty=2)
		if(focal_row) points(D[focal_row], CD[focal_row], col='red')
	
		plot(CA ~ CD, xlab='Crown diameter [m]', ylab='Crown projection area [m2]', type='l')
		if(focal_row) points(Cd[focal_row], CA[focal_row], col='red')
	
		plot(LAI ~ D, xlab='Stem diameter [m]', ylab='Leaf area index', type='l')
		abline(v=Dmax, col='grey', lty=2)
		if(focal_row) points(D[focal_row], LAI[focal_row], col='red')

		plot(BT ~ D, xlab='Stem diameter [m]', ylab='Biomass [t_ODM]', type='l')
		abline(v=Dmax, col='grey', lty=2)
		if(focal_row) points(D[focal_row], BT[focal_row], col='red')
		
	})
}

plot_tree <- function(data, age=NULL){
	
	# visualise a tree at a given age
	
	# - get the plot window for the tree
	CD_max = max(data$CD)
	H_max = max(data$H)	
	if(is.null(age)) row <- max(data$AGE)
	plot.new()
	plot.window(c(-CD_max, CD_max), c(0, H_max),  asp = 1.0)
	box()
	axis(1)
	axis(4)
	
	# Simple boxes to represent the trunk and canopy
	with(data[row,],{
		rect(-CD/2, H - (H*CLP), CD/2, H, col='#005500CC')
		rect(-D/2, 0, D/2, H, col='grey')
	})
}

plot_geometry_overlay <- function(res, settings){
	
	# Takes the simulation data from a formind res file
	# for a single tree and overplots geometry predictions 
	# for given settings using the R functions
	
	plot(H ~ D, data=res, type='l', lwd=3, col='grey')
	res$HR <- D_to_H(res$D, settings) 
	lines(HR ~ D, data=res, col='red')
	deltaH <- with(res, HR - H)
	plot(density(deltaH))
	
	plot(CLP ~ H, data=res, type='l', lwd=3, col='grey')
	res$CLPR <- H_to_Clp(res$H, settings) 
	lines(CLPR ~ H, data=res, col='red')
	deltaClp <- with(res, CLPR - CLP)
	plot(density(deltaClp))
	
	plot(CD ~ D, data=res, type='l', lwd=3, col='grey')
	res$CDR <- D_to_Cd(res$D, settings)
	lines(CDR ~ D, data=res, col='red')
	deltaCd <- with(res, CDR - CD)
	plot(density(deltaCd))

	plot(LAI ~ D, data=res, type='l', lwd=3, col='grey')
	res$LAIR <- D_to_LAI(res$D, settings)
	lines(LAIR ~ D, data=res, col='red')
	deltaLAI <- with(res, LAIR - LAI)
	plot(density(deltaLAI))

	plot(BT ~ D, data=res, type='l', lwd=3, col='grey')
	res$BTR <- Geom_to_B(res$D, res$H, settings)
	lines(BTR ~ D, data=res, col='red')
	deltaBT <- with(res, BTR - BT)
	plot(density(deltaBT))

}

plot_simulation_overlay <- function(res, sim, vars=c('D', 'H', 'CD', 'CLP', 'BT' ,
									'LAI', 'PB', 'DInc', 'BInc', 'R_Growth', 'R_Main'), ...){

	# Takes a res file and a simulation data frame run over
	# the same time frame and plots the differences between
	# the variables by Age to show mismatch
	
	for(v in vars){
		plot((res[[v]] - sim[[v]]) ~ res$AGE, type='l', main=v, ylab=v, xlab='Age', ...)
	}	
	
}