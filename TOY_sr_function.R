#------------------------------------------------------------------------------#
# Spawner-recruitment simulation function
#------------------------------------------------------------------------------#
# ny <- number of brood years
# alpha <- intrinsic productivity (not in log space)
# beta <- density dependence 
# sigma.R <- log-normal recruitment variation
# phi <- lag-one correlation in recruitment through time
# mat <- maturation schedules (assumes max age is 7 years and fish mature at 4-7 years) ##vector length 4
# U_alpha <- alpha parameter for beta distribution for annual harvest rate (see https://keisan.casio.com/exec/system/1180573226)
# U_beta <- beta parameter for beta distribution for annual harvest rate

#read in values to dissect components of the function
if(FALSE){
  sim_sr <- process(ny = 40,
                    alpha = 4,
                    beta = 0.0005,
                    phi = 0.6,
                    mat = c(0.2, 0.3, 0.4, 0.1),
                    sigma.R = 0.5,
                    U_alpha = 7,
                    U_beta = 7)
}


process = function(ny, phi, mat, alpha, beta, sigma.R, U_alpha, U_beta){
  ns = length(alpha) # number of stocks
  epi = rnorm(ny, sd = sigma.R) # recruitment deviations
  Ro = log(alpha)/beta # starting spawner abundances
  R = t(matrix(0,ns,ny)) # placeholder for recruitment time-series 
  S = R # placeholder for spawner time-series 
  v = R # placeholder for component of recruitment deviations that is correlated through time 
  N = array(0,dim=c(ny,4,ns)) # placeholder for returns by age (4-7 yrs) 
  Ntot = R # placeholder for total returns
  H = Ntot # placeholder for harvest
  S = Ntot # placeholder for spawners
  predR = Ntot # placeholder for predicted recruitment
  
  R[1:7,] = t(replicate(7, Ro, simplify=T)) * exp(epi[1:7]) # recruitment in first 7 years without corresponding brood year spawning abundance ##why transpose
  
  # Loop through years of simulation	  
  for(i in (7+1):ny){ 
    
    # returns in year i are recruits of age X from brood year X years prior
    N[i,1,1] = R[i-(4),1] * mat[1]
    N[i,2,1] = R[i-(5),1] * mat[2]
    N[i,3,1] = R[i-(6),1] * mat[3]
    N[i,4,1] = R[i-(7),1] * mat[4]
    
    Ntot[i,1] = sum(N[i,,1])
    
    # apply harvest 
    H[i,1] =  rbeta(1,U_alpha, U_beta)*Ntot[i,1]
    S_exp = Ntot[i,1]-H[i,1] ; S_exp[S_exp<0] = 0
    S[i,1] = S_exp
    
    # predict recruitment
    R[i,1] = alpha*S[i,1]*exp(-beta*S[i,1]+phi*v[i-1,1]+epi[i])
    predR[i,1] = alpha*S[i,1]*exp(-beta*S[i,1])
    v[i,1] = log(R[i,1])-log(predR[i,1])
    v[v[,1]=='NaN'] <- 0
  }
  
  list(S=S[,], R=R[,], N=Ntot[,], N.age=N[,,1]) ##would be clearer if the N.age notation wasn't based on position (e.g. 1,2,3,4) but age (e.g. 4,5,6,7)
}
