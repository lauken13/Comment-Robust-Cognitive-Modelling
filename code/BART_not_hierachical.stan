data {
  int<lower=1> nconds;
  int<lower=1> ntrials;
  vector<lower=0,upper=1>[3] p;
  int<lower=1,upper=3> p_balloon[ntrials]; //p value for each trial. 
  int<lower=1> options[nconds,ntrials];
  int d[nconds,ntrials,30];
  int<lower=0> up_prior;
}
parameters {
  real<lower=0, upper=up_prior> gplus;
  real<lower=0, upper=up_prior> beta;
} 
transformed parameters {
  vector<lower=0>[nconds] gplus_rep = rep_vector(gplus,3);
  vector<lower=0>[nconds] beta_rep = rep_vector(beta,3);
  matrix<lower=0>[nconds,3] omega;
  // Optimal Number of Pumps
  for (i in 1:3){
      omega[,i] = -gplus_rep / log1m(p[i]);
  }
}
model {
  target += uniform_lpdf(gplus | 0,up_prior);
  target += uniform_lpdf(beta | 0,up_prior);
//  target += normal_lpdf(gplus_cons | mug, sigmag) - normal_lccdf(0 | mug, sigmag);
//  target += normal_lpdf(beta_cons | mub, sigmab) - normal_lccdf(0 | mub, sigmab);
  // Choice Data
  for (i in 1:nconds) {
    for (j in 1:ntrials) {
      for (k in 1:options[i,j]) {
        target += bernoulli_logit_lpmf(abs(d[i,j,k]-1) | -beta_rep[i] * (k - omega[i,p_balloon[j]]));
      }
    }
  }
}
generated quantities {
  real log_lik[nconds,ntrials,30];
  for (i in 1:nconds) {
    for (j in 1:ntrials) {
      for (k in 1:options[i,j]) {
        log_lik[i,j,k] = bernoulli_logit_lpmf(abs(d[i,j,k]-1) | -beta_rep[i] * (k - omega[i,p_balloon[j]]));
      }
    }
  }
}
