source('../FORMIND_tree_models.R')

## GEOMETRY Check
settings <- list(Hmax = 35.0,
				 D_to_H = list(method='power', pars=list(h_0=4.5 * 100^0.45, h_1=0.45)),
				 D_to_Cd = list(method='power', pars=list(cd_0=0.55*100^0.7, cd_1=0.7, cd_2=0.15)),
				 H_to_Clp = list(method='linear', pars=list(cl_0=0.35)),
				 Geom_to_B = list(method='geometric', pars=list(rho=0.4, sigma=0.7, f_0=0.775, f_1=-0.175),
								  form_factor=16),
				 D_to_LAI = list(method='power', pars=list(l_0=2, l_1=0.2)),
				 D_to_Dinc = list(method='polynomial', pars=list(dmin=0, gmax=0.04, dpropgmax=1/3, gstart=0.4, gend=0.1)),
				 GPP = list(pars=list(alpha=0.15, pmax=15.0, k=0.6, m=0.1, active=365)),
				 r_g = 0.2)

settings <- set_maxima(settings)
curves <- geometry_curves(settings)

# look at the curves generated in R
par(mfrow=c(3,2), mar=c(3,3,1,1), mgp=c(2,0.8,0))
plot_geometry(curves)

# Comparision to output of lonesome.par, which is supposed to contain a single instance of 
# a tree growing with these settings.
tree <- read.delim('lonesome.res', skip=2)
# drop data after the tree reaches maximum size
tree <- tree[1:which.max(tree$H),]

# Check the Geometry by overlaying R geom on res outputs
par(mfrow=c(5,2), mar=c(3,3,1,1), mgp=c(2,0.8,0))
plot_geometry_overlay(tree, settings)

## PRODUCTION
production <- growth_simulation(settings, Dinit=0.035, n_years=56)
production2 <- growth_simulation(settings, Dinit=0.035, n_years=56, truncate_overgrowth=TRUE)


pdf('Lonesome_growth.pdf', paper='a4', height=8, width=8)
	par(mfrow=c(2,2), mar=c(3,3,1,1), mgp=c(2,0.8,0))
	plot_production(tree, main='FORMIND .res file contents')
	plot_production(production, ylim=c(-0.1,1.3), main='R (not constrained)')
	plot_production(production2, main='R (constrained)')

	# This one isn't the same because the polynomial increment calculation is off
	plot(DInc ~ D, data=production, type='l', col='red')
	lines(DInc ~ D, data=tree, type='l', lwd=3, col='grey')
	lines(DInc ~ D, data=production2, type='l', col='red', lty=2)
	legend('topright', lty=c(1,1,2), lwd=c(3,1,1), col=c('grey','red','red'), 
			legend=c('FORMIND .res', 'R (unconstrained)', 'R (constrained)'),
			bty='n')
dev.off()