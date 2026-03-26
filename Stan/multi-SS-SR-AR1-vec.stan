// Stan version of age-structured state-space spawner-recruitment model with n_aR-1 process variation (adapted from Fleischman et al. CJFn_aS. 2013)

data{
  int n_pops;             // number of populaitons
  int n_Y;                // number of calendar years of obs for all pops
  int n_Y_R;              // number of recruitment years for all pops 
  int n_y;                // single pop obs length
  int n_y_R;              // single pop rec length 
  vector[n_Y] ind_Y;      // vector indexing data to populations
  vector[n_Y_R] ind_Y_R;  // vector indexing recruitment estimates to populations
  int a_min;              // minimum age class
  int a_max;              // maximum age class
  int n_a;                // number of age classes
  int A_obs[n_y, n_a];    // observed age composition in counts by age class
  vector[n_Y] S_obs;      // observed spawners
  vector[n_Y] H_obs;      // observed harvest
  vector[n_Y] S_CV;       // spawner observation error CV
  vector[n_Y] H_CV;       // harvest observation error CV
}

parameters{
  vector<lower=0>[n_Y_R] lnR;             // log recruitment states
  vector<lower=0>[n_pops] lnalpha;        // log Ricker a
  vector<lower=0>[n_pops] beta;           // Ricker b
  vector<lower=0>[n_pops]sigma_R;         // process error
  vector<lower=0>[n_pops] sigma_R0;       // process error for first a.max years with no spawner link
  vector<lower=-1, upper=1>[n_pops] phi;  // lag-1 correlation in process error
  vector[n_pops] lnresid_0;               // first residual for lag-1 correlation in process error
  real<lower=0> mean_ln_R0;               // "true" mean log recruitment in first a.max years with no spawner link
  vector<lower=0.01, upper=0.99>[n_Y] U;  // harvest rate
  vector<lower=0, upper=1>[n_a - 1] prob; // maturity schedule probs
  real<lower=0, upper=1> D_scale;         // governs variability of age proportion vectors across cohorts
  matrix<lower=0.01>[n_y_R, n_a] g;       // individual year/age class gamma variates for generating age at maturity proportions
}

transformed parameters{
  vector<lower=0>[n_Y] N;                 // run size states
  vector<lower=0>[n_Y] S;                 // spawner states
  vector[n_Y] lnS;                        // log spawner states
  vector<lower=0>[n_Y] C;                 // catch states
  vector[n_Y] lnC;                        // log catch states
  vector<lower=0>[n_Y_R] R;               // recruitment states
  vector[n_Y_R] lnresid;                  // log recruitment residuals
  vector[n_Y_R] lnRm_1;                   // log recruitment states in absence of lag-one correlation
  vector[n_Y_R] lnRm_2;                   // log recruitment states after accounting for lag-one correlation
  matrix<lower=0>[n_Y, n_a] N_ta;         // returns by age matrix
  matrix<lower=0, upper=1>[n_y_R, n_a] p; // age at maturity proportions
  vector<lower=0,upper=1>[n_a] pi;        // maturity schedule probs
  real<lower=0> D_sum;                    // inverse of D_scale which governs variability of age proportion vectors across cohorts
  vector<lower=0>[n_a] Dir_alpha;         // Dirichlet shape parameter for gamma distribution used to generate vector of age-at-maturity proportions
  matrix<lower=0, upper=1>[n_y, n_a] q;   // age composition by year/age classr matrix
  vector<lower=0>[n_pops] S_max;          // Spawners that maximize recruitment  

  // Maturity schedule: use a common maturation schedule to draw the brood year specific schedules
  pi[1] = prob[1];
  pi[2] = prob[2] * (1 - pi[1]);
  pi[3] = prob[3] * (1 - pi[1] - pi[2]);
  pi[4] = 1 - pi[1] - pi[2] - pi[3];
  D_sum = 1/D_scale^2;

  for(a in 1:n_a){
    Dir_alpha[a] = D_sum * pi[a];
    for(y in 1:n_y_R){
      p[y, a] = g[y, a]/sum(g[y, ]);
    }
  }

  R = exp(lnR);
  lnS[t] = log(S[t]);    
  lnC[t] = log(C[t]);
  S_max = 1/beta;

  // Calculate the numbers at age matrix as brood year recruits at age (proportion that matured that year)
  // THIS is the tricky part, we want to keep the proportions the same for pops, 
    // but they need to be partitioned into their specific N_ta
    //I guess the best way is to figure out how to index Nta and R here to be long
  //for (t in 1:n_Y) {
  //  for(a in 1:n_a){
  //    N_ta[t, a] = R[t + n_a - a] * p[t + n_a - a, a]; // p has to stay n_y but we need bigger R and N_ta 
  //  }
  //}
  for(p in 1:n_pops){ //1:2
    for(t in 1:n_y){ //1:33
      for(a in 1:n_a){ //1:4
        i = n_y*(p-1)+t //new i index to make N_ta and R long; so stupid it might work... 
        N_ta[i, a] = R[i + n_a - a] * p[t + n_a - a, a];// but this doesnt help b/c p is wrong dim still
      }
    }
  }

  // Calculate returns, spawners and catch by return year
    // does this need to be a loop?? it's all vectors... 
  for(t in 1:n_Y){
    N[t] = sum(N_ta[t, 1:n_a]);
    S[t] = N[t] * (1 - U[t]);
    C[t] = N[t] * U[t];
  }

  // Calculate age proportions by return year
  for(t in 1:n_y){
    for(a in 1:n_a){
      q[t,a] = N_ta[t, a]/N[t]; //BAD - can't make bigger obj smaller into q
    }
  }

  // Ricker SR with AR1 process on log recruitment residuals for years with brood year spawners
  for (i in 1:n_y_R) {
    lnresid[i] = 0.0;
    lnRm_1[i] = 0.0;
    lnRm_2[i] = 0.0;
  }

  for (y in (n_a + a_min):n_y_R) {
    lnRm_1[y] = lnS[y - a_max] + lnalpha - beta * S[y - a_max];
    lnresid[y] = lnR[y] - lnRm_1[y];
  }

  lnRm_2[n_a + a_min] =  lnRm_1[n_a + a_min] + phi * lnresid_0;

  for (y in (n_a + a_min + 1):n_y_R) {
    lnRm_2[y] =  lnRm_1[y] + phi * lnresid[y - 1];
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
  for (y in 1:n_y_R) {
    for (a in 1:n_a) {
      //g[y,a] ~ gamma(Dir_alpha[a],1);
      target += gamma_lpdf(g[y, a]|Dir_alpha[a], 1);
    }
  }

  // First `a.max` years of recruits, for which there is no spawner link
  lnR[1:a_max] ~ normal(mean_ln_R0, sigma_R0);

  // State model
  lnR[(n_a + a_min):n_y_R] ~ normal(lnRm_2[(n_a + a_min):n_y_R], sigma_R);

  // Observation model
  for(t in 1:n_y){
    target += multinomial_lpmf(A_obs[t,1:n_a]|to_vector(q[t,1:n_a]));
    }
    
  for(t in 1:n_Y){
    U[t] ~ beta(1,1);
    H_obs[t] ~ lognormal(lnC[t], sqrt(log((H_CV[t]^2) + 1)));
    S_obs[t] ~ lognormal(lnS[t], sqrt(log((S_CV[t]^2) + 1)));
  }
}
