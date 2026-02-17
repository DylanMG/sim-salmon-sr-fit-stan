#------------------------------------------------------------------------------#
# Spawner-recruitment simulation function
#------------------------------------------------------------------------------#
#n_y <- number of brood years
#alpha <- intrinsic productivity (not in log space)
#beta <- density dependence 
#sigma_R <- log-normal recruitment variation
#phi <- lag-one correlation in recruitment through time
#mat <- maturation schedules (assumes max age is 7 years and fish mature at 4-7 years)
#U_alpha <- alpha parameter (i.e. shape) for beta distribution for annual harvest rate (see 
#U_beta <- beta parameter ' '

process = function(n_y, a_max, phi, mat, alpha, beta, sigma_R, U_alpha, U_beta){
  n_s <- length(alpha) #number of stocks
  eps <- rnorm(n_y, sd = sigma_R) #epsilon - recruitment deviations
  R0 <- log(alpha)/beta #unobserved, starting spawner abundances (i.e. at equilibrium)
  R <- t(matrix(0, n_s, n_y)) #placeholder for recruitment time-series 
  S <- R #placeholder for spawner time-series 
  v <- R #placeholder for component of recruitment deviations that is correlated through time 
  N <- array(0, dim = c(n_y, 4, n_s)) #placeholder for returns by age (4-7 yrs) 
  N_tot <- R #placeholder for total returns
  H <- N_tot #placeholder for harvest
  S <- N_tot #placeholder for spawners
  R_pred <- N_tot #placeholder for predicted recruitment
  
  R[1:a_max,] <- t(rep(R0, a_max)) * exp(eps[1:a_max]) #recruits with noise in first not fully observed years
  
  # Loop through "fully observed" years of simulation	  
  for(i in (a_max+1):n_y){ 
    
    # returns in year i are recruits of age X from brood year X years prior
    N[i, 1, 1] <- R[i - (4), 1] * mat[1]
    N[i, 2, 1] <- R[i - (5), 1] * mat[2]
    N[i, 3, 1] <- R[i - (6), 1] * mat[3]
    N[i, 4, 1] <- R[i - (7), 1] * mat[4]
    
    N_tot[i, 1] <- sum(N[i, , 1])
    
    # apply harvest 
    H[i, 1] <-  rbeta(1, U_alpha, U_beta) * N_tot[i, 1]
    S_exp <- N_tot[i, 1] - H[i, 1] 
    S_exp[S_exp < 0] <- 0
    S[i, 1] <- S_exp
    
    # predict recruitment
    R[i, 1] <- alpha * S[i, 1] * exp(-beta * S[i, 1] + phi * v[i - 1, 1]+eps[i])
    R_pred[i, 1] <- alpha* S[i, 1] * exp(-beta * S[i, 1])
    v[i, 1] <- log(R[i, 1]) - log(R_pred[i, 1])
    v[v[, 1] == 'NaN'] <- 0
  }
  
  list(S = round(S[ , ], 0), N_age = round(N[ , , 1], 0), N = round(N_tot[ , ], 0), 
       R = round(R[ , ], 0))
}
