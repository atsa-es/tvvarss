data {
  int<lower=0> n_year;
  int<lower=0> n_site;
  int<lower=0> n_spp;
  int<lower=0> b_diag[n_spp*n_spp]; # vector indicating whether elements are on diagonal
  matrix[n_site,n_spp] x0; # data
  int<lower=0> shared_q[n_spp,n_site+1]; # matrix indicating which sites/spp share residual process variances
  int<lower=0> n_q; # number of variances being estimated, max(shared_q)
  int<lower=0> shared_r[n_spp,n_site+1]; # matrix indicating which sites/spp share obs variances
  int<lower=0> n_r; # number of obs variances being estimated, max(shared_q)
  int<lower=0> shared_u[n_spp,n_site+1]; # matrix indicating which sites/spp share trends
  int<lower=0> n_u; # number of trends being estimated, max(shared_u)
  int<lower=0> est_trend;
  int<lower=0> demean;
  int<lower=0> row_indices[(n_spp*n_spp)];
  int<lower=0> col_indices[(n_spp*n_spp)];
  int<lower=0> n_pos;
  int<lower=0> spp_indices_pos[n_pos];
  int<lower=0> site_indices_pos[n_pos];
  int<lower=0> year_indices_pos[n_pos];
  real y[n_pos]; # vector[n_pos] y; # data
}
parameters {
  vector[(n_spp*n_spp)] vecB[n_year]; # elements accessed [n_year,n_spp]
  real<lower=0> sigma_rw_pars[2]; # sds for random walk
  matrix[n_year,n_spp] x[n_site]; # unobserved states
  real<lower=0> resid_process_sd[n_q]; # residual sds
  real<lower=0> obs_sd[n_r]; # residual sds
  real u[n_u]; # trends
}
transformed parameters {
  vector<lower=0>[(n_spp*n_spp)] sigma_rw;
  matrix<lower=0>[n_spp, n_site] resid_process_mat;
  matrix<lower=0>[n_spp, n_site] obs_mat;
  matrix[n_spp, n_site] u_mat;
  matrix[n_spp,n_spp] B[(n_year-1)]; # B matrix, accessed as n_year, n_spp, n_spp
  matrix[n_year,n_spp] pred[n_site]; # predicted unobserved states

  for(i in 1:(n_spp*n_spp)) {
    sigma_rw[i] = sigma_rw_pars[b_diag[i]];
  }
  for(i in 1:n_spp) {
    for(j in 1:n_site) {
      resid_process_mat[i,j] = resid_process_sd[shared_q[i,j]];
      u_mat[i,j] = u[shared_u[i,j]];
      obs_mat[i,j] = obs_sd[shared_r[i,j]];
    }
  }

  for(s in 1:n_site) {
    pred[s,1,] = x[s,1,]; # states for first year
  }
  for(t in 2:n_year) {
    # fill in B matrix, shared across sites
    for(i in 1:(n_spp*n_spp)) {
      B[t-1,row_indices[i],col_indices[i]] = vecB[t-1,i];
    }

    # do projection to calculate predicted values. modify code depending on whether
    # predictions should be demeaned before projected, and whether or not trend included.
    for(s in 1:n_site) {
     if(est_trend == 0) {
      if(demean==0) {pred[s,t,] = x[s,t-1,] * B[t-1,,];}
      if(demean==1) {pred[s,t,] = (x[s,t-1,]-pred[s,t-1,]) * B[t-1,,];}
     }
     if(est_trend == 1) {
      if(demean==0) {pred[s,t,] = x[s,t-1,] * B[t-1,,] + u_mat[,s]';}
      if(demean==1) {pred[s,t,] = (x[s,t-1,]-pred[s,t-1,]) * B[t-1,,] + u_mat[,s]';}
     }
    }
  }
}
model {
  sigma_rw_pars[1] ~ cauchy(0,5);
  sigma_rw_pars[2] ~ cauchy(0,5);
  vecB[1] ~ normal(0, 3); # prior for first state
  for(t in 2:n_year) {
    vecB[t] ~ normal(vecB[t-1], sigma_rw); # vectorized random in B
  }
  # prior on first time step
  for(site in 1:n_site) {
    for(spp in 1:n_spp) {
      x[site,1,spp] ~ normal(x0[site,spp],1);
    }
  }
  # process model for remaining sites
  for(t in 2:n_year) {
  for(site in 1:n_site) {
    for(spp in 1:n_spp) {
      x[site,t,spp] ~ normal(pred[site,t,spp], resid_process_mat[spp,site]);
    }
  }
  }
  for(i in 1:n_q) {
    # prior on process standard deviations
    resid_process_sd[i] ~ cauchy(0,5);
  }
  for(i in 1:n_u) {
    # prior on trends
    u[i] ~ normal(0,1);
  }

  # gaussian likelihood for now
  for(i in 1:n_pos) {
    y[i] ~ normal(x[site_indices_pos[i],year_indices_pos[i],spp_indices_pos[i]], obs_mat[spp_indices_pos[i],site_indices_pos[i]]);
  }

}
generated quantities {
  vector[n_pos] log_lik;
  # for use in loo() package
  for (n in 1:n_pos) log_lik[n] = normal_lpdf(y[n] | x[site_indices_pos[n],year_indices_pos[n],spp_indices_pos[n]], obs_mat[spp_indices_pos[n],site_indices_pos[n]]);
}
