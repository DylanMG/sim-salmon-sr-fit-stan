// Stan version of age-structured state-space spawner-recruitment model with n_aR-1 process variation (adapted from Fleischman et al. CJFn_aS. 2013)

data{
  int n_y;                  // number of calender years
  int n_pops;               // number of populations
  int a_min;                // minimum age class
  int a_max;                // maximum age class
  int n_a;                  // number of age classes
  int n_y_R;                // number of recruitment years
  int A_obs[n_y, n_a];      // observed age composition in counts by age class
  matrix[n_y, n_pops] S_obs;// observed spawners
  matrix[n_y, n_pops] H_obs;// observed harvest
  vector[n_y] S_CV;         // spawner observation error CV
  vector[n_y] H_CV;         // harvest observation error CV
}

parameters{
  matrix<lower=0>[n_y_R, n_pops] lnR;     // log recruitment states
  vector<lower=0>[n_pops] lnalpha;        // log Ricker a
  vector<lower=0>[n_pops] beta;           // Ricker b
  vector<lower=0>[n_pops]sigma_R;         // process error
  vector<lower=0>[n_pops] sigma_R0;       // process error for first a.max years with no spawner link
  vector<lower=-1, upper=1>[n_pops] phi;  // lag-1 correlation in process error
  vector[n_pops] lnresid_0;               // first residual for lag-1 correlation in process error
  vector<lower=0>[n_pops] mean_ln_R0;     // "true" mean log recruitment in first a.max years with no spawner link
  matrix<lower=0.01, upper=0.99>[n_y, n_pops] U;// harvest rate
  //shared
  vector<lower=0, upper=1>[n_a - 1] prob; // maturity schedule probs
  real<lower=0, upper=1> D_scale;         // governs variability of age proportion vectors across cohorts
  matrix<lower=0.01>[n_y_R, n_a] g;       // individual year/age class gamma variates for generating age at maturity proportions
}

transformed parameters{
  matrix<lower=0>[n_y, n_pops] N;         // run size states
  matrix<lower=0>[n_y, n_pops] S;         // spawner states
  matrix[n_y, n_pops] lnS;                // log spawner states
  matrix<lower=0>[n_y, n_pops] C;         // catch states
  matrix[n_y, n_pops] lnC;                // log catch states
  matrix<lower=0>[n_y_R, n_pops] R;       // recruitment states
  matrix[n_y_R, n_pops] lnresid;          // log recruitment residuals
  matrix[n_y_R, n_pops] lnRm_1;           // log recruitment states in absence of lag-one correlation
  matrix[n_y_R, n_pops] lnRm_2;           // log recruitment states after accounting for lag-one correlation
  real<lower=0>N_ta[n_y, n_a, n_pops];    // returns by age matrix
  //shared
  matrix<lower=0, upper=1>[n_y_R, n_a] p; // age at maturity proportions
  vector<lower=0,upper=1>[n_a] pi;        // maturity schedule probs
  real<lower=0> D_sum;                    // inverse of D_scale which governs variability of age proportion vectors across cohorts
  vector<lower=0>[n_a] Dir_alpha;         // Dirichlet shape parameter for gamma distribution used to generate vector of age-at-maturity proportions
  matrix<lower=0, upper=1>[n_y, n_a] q;   // age composition by year/age class matrix
  vector<lower=0>[n_pops] S_max;          // Spawners that maximize recruitment  

  // Maturity schedule: use a common maturation schedule to draw the brood year specific schedules
  pi[1] = prob[1];
  pi[2] = prob[2] * (1 - pi[1]);
  pi[3] = prob[3] * (1 - pi[1] - pi[2]);
  pi[4] = 1 - pi[1] - pi[2] - pi[3];
  D_sum = 1/D_scale^2;

  for(a in 1:n_a){
    Dir_alpha[a] = D_sum * pi[a];
    for (y in 1:n_y_R) {
      p[y, a] = g[y, a]/sum(g[y, ]);
    }
  }

  R = exp(lnR);
  S_max = 1/beta;

  // Calculate the numbers at age matrix as brood year recruits at age (proportion that matured that year)
  for(i in 1:n_pops){
    for(t in 1:n_y){
      for(a in 1:n_a){
        N_ta[t, a, i] = R[t + n_a - a, i] * p[t + n_a - a, a];
    }
   }
  }

  // Calculate returns, spawners and catch by return year
  for(i in 1:n_pops){
    for(t in 1:n_y){
    N[t, i] = sum(N_ta[t, 1:n_a, i]);
    S[t, i] = N[t, i] * (1 - U[t, i]);
    lnS[t, i] = log(S[t, i]);
    C[t, i] = N[t, i] * U[t, i];
    lnC[t, i] = log(C[t, i]);
    }
  }
  
  // Calculate age proportions by return year
    for(t in 1:n_y){
      for(a in 1:n_a){
        q[t,a] = sum(N_ta[t, a, ])/sum(N[t, ]); //THIS is where the multi-stock dynamics get "compressed" into 1 q
        }
      }

  // Ricker SR with AR1 process on log recruitment residuals for years with brood year spawners
  for(i in 1:n_pops){
    for(j in 1:n_y_R){
    lnresid[j, i] = 0.0;
    lnRm_1[j, i] = 0.0;
    lnRm_2[j, i] = 0.0;
    }
  }
  
  for(i in 1:n_pops){
    for(y in (n_a + a_min):n_y_R){
    lnRm_1[y, i] = lnS[y - a_max, i] + lnalpha[i] - beta[i] * S[y - a_max, i];
    lnresid[y, i] = lnR[y, i] - lnRm_1[y, i];
    }
    lnRm_2[n_a + a_min, i] =  lnRm_1[n_a + a_min, i] + phi[i] * lnresid_0[i];
  }
  
  for(i in 1:n_pops){
    for(y in (n_a + a_min + 1):n_y_R){
    lnRm_2[y, i] =  lnRm_1[y, i] + phi[i] * lnresid[y - 1, i];
    }
  }
}

model{
  // Priors
  lnalpha ~ normal(0,3);
  beta ~ normal(0,1);
  sigma_R ~ normal(0,2);
  lnresid_0 ~ normal(0,20);
  mean_ln_R0 ~ normal(0,20);
  sigma_R0 ~ inv_gamma(2,1); 
  prob[1] ~ beta(1,1);
  prob[2] ~ beta(1,1);
  prob[3] ~ beta(1,1);
  D_scale ~ beta(1,1);

  // Likelihoods
  // Gamma variates for each year and age class which are used to determine age at maturity proportions
  for(y in 1:n_y_R){
    for(a in 1:n_a){
      //g[y,a] ~ gamma(Dir_alpha[a],1);
      target += gamma_lpdf(g[y, a]|Dir_alpha[a], 1);
    }
  }
  
  for(i in 1:n_pops){
    // First `a.max` years of recruits, for which there is no spawner link
    lnR[1:a_max, i] ~ normal(mean_ln_R0[i], sigma_R0[i]);

    // State model
    lnR[(n_a + a_min):n_y_R, i] ~ normal(lnRm_2[(n_a + a_min):n_y_R, i], sigma_R[i]);
  }
  
  // Observation model
  for(t in 1:n_y){
  //A_obs[t,1:n_a]) ~ multinomial(q[t,1:n_a]);
    target += multinomial_lpmf(A_obs[t,1:n_a]|to_vector(q[t,1:n_a]));
    for(i in 1:n_pops){
      U[t, i] ~ beta(1,1);
      H_obs[t, i] ~ lognormal(lnC[t, i], sqrt(log((H_CV[t]^2) + 1)));
      S_obs[t, i] ~ lognormal(lnS[t, i], sqrt(log((S_CV[t]^2) + 1)));
    }
  }
}
