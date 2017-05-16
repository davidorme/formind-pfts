## Settings object - example settings using the tropical forest project parameters
## for three functional types

source('FORMIND_tree_models.R')

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

# look at the geometry and growth
pioneer <- growth_simulation(PFTs[['pioneer']], Dinit=0.001, n_years=145)
par(mfrow=c(3,2), mar=c(3,3,1,1), mgp=c(2,0.8,0))
plot_geometry(pioneer)
plot_production(pioneer)

intermed <- growth_simulation(PFTs[['intermed']], Dinit=0.001, n_years=817)
par(mfrow=c(3,2), mar=c(3,3,1,1), mgp=c(2,0.8,0))
plot_geometry(intermed)
plot_production(intermed)

climax <- growth_simulation(PFTs[['climax']], Dinit=0.001, n_years=2001)
par(mfrow=c(3,2), mar=c(3,3,1,1), mgp=c(2,0.8,0))
plot_geometry(climax)
plot_production(climax)

