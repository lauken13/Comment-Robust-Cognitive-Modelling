data {
  int<lower=1> nconds;
  int<lower=1> ntrials;
  vector<lower=0,upper=1>[3] p; // Three different p values
  int<lower=1,upper=3> p_balloon[ntrials]; //p value for each trial. 
  int<lower=1> options[nconds,ntrials];
  int d[nconds,ntrials,30];
  int<lower=0> up_prior;
}
parameters {
  vector<lower=0>[nconds] gplus_raw;
  vector<lower=0>[nconds] beta_raw;
  // Priors
  real<lower=0,upper=up_prior> mug;
  real<lower=0,upper=up_prior> sigmag;
  real<lower=0,upper=up_prior> mub;
  real<lower=0,upper=up_prior> sigmab;
} 
transformed parameters {
  matrix<lower=0>[nconds,3] omega;
  vector<lower=0>[nconds] gplus;
  vector<lower=0>[nconds] beta;
  // Optimal Number of Pumps
  gplus = mug + gplus_raw*sigmag;
  beta = mub + beta_raw*sigmab;
  for (i in 1:3){
      omega[,i] = -gplus / log1m(p[i]);
  }
}
model {
  target += uniform_lpdf(mug | 0, up_prior);
  target += uniform_lpdf(sigmag | 0,up_prior);
  target += uniform_lpdf(mub | 0,up_prior);
  target += uniform_lpdf(sigmab | 0,up_prior);
  target += normal_lpdf(gplus_raw | 0,1) - normal_lccdf(0|0,1);
  target += normal_lpdf(beta_raw | 0,1) - normal_lccdf(0|0,1);
//  for (i in 1:nconds) {
//    target += normal_lpdf(gplus[i] | mug, 1/sigmag) - normal_lccdf(0 | mug, sigmag);
//    target += normal_lpdf(beta[i] | mub, sigmab) - normal_lccdf(0 | mub, sigmab);
//    }
  // Choice Data
  for (i in 1:nconds) {
    for (j in 1:ntrials) {
      for (k in 1:options[i,j]) {
        target += bernoulli_logit_lpmf(abs(d[i,j,k]-1) | -beta[i] * (k - omega[i,p_balloon[j]]));
      }
    }
  }
}
generated quantities {
  real log_lik[nconds,ntrials,30];
  for (i in 1:nconds) {
    for (j in 1:ntrials) {
      for (k in 1:options[i,j]) {
        log_lik[i,j,k] = bernoulli_logit_lpmf(abs(d[i,j,k]-1) | -beta[i] * (k - omega[i,p_balloon[j]]));
      }
    }
  }
}

