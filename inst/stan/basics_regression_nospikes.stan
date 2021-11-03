functions {
  matrix designMatrix(int l, 
                      vector mu, 
                      real rbf_variance,
                      vector ml,
                      int q) {

    // vector [q] x = log(mu);
    vector [q] x = mu;
    real h;
    matrix [q, l] X;
    h = (ml[2] - ml[1]) * rbf_variance;
    for (i in 1:q) {
      for (j in 1:l) {
        X[i, j] = 1;
      }
    }
    X[, 2] = x;
    for (i in 1:(l - 2)) {
      vector[q] tmp;
      for (j in 1:q) {
        tmp[j] = pow(x[j] - ml[i], 2);
      }
      X[, i + 2] = exp(-0.5 * tmp / pow(h, 2));
    }
    return(X);
  }
}

data {
  int q;
  int n;
  int p;
  int<lower=0> counts[q, n];
  real as;
  real bs;
  real atheta;
  real btheta;
  matrix [n, p] batch_design;
  vector [n] aphi;
  vector[q] mu_mu;
  real smu;
  real mu0; // fixed mean log expression
  real astwo;
  real bstwo;
  int l;
  int refgene;
  int notrefgene [q - 1];
  vector [n] size_factors;
  vector [l] mbeta;
  matrix[l, l] vbeta;
  real rbf_variance;
  vector[l - 2] ml;
  real eta;
}

transformed data {
  vector [q-1] mu_mean = rep_vector(mu0, q - 1);
  vector [q-1] q_ones = rep_vector(1, q - 1);
  cov_matrix [q-1] mu_cov = smu * (diag_matrix(q_ones) - (q_ones * transpose(q_ones)) / q);
}


parameters {
  vector [q] mu;
  vector <lower=0> [q] delta;
  vector [l] beta;
  real <lower=0> stwo;
  vector <lower=0> [q] lambda;
}

transformed parameters {
  vector [q] fu = designMatrix(l, mu, rbf_variance, ml, q) * beta;
  vector [q] epsilon = log(delta) - fu;
  real sum_mu = sum(mu[notrefgene]);
}

model {
  stwo ~ inv_gamma(astwo, bstwo);
  beta ~ multi_normal(mbeta, stwo * vbeta);
  mu[notrefgene] ~ multi_normal(mu_mean, mu_cov);
  mu[refgene] ~ normal((q * mu0) - sum_mu, 0.01);
  lambda ~ gamma(eta / 2, eta / 2);
  delta ~ lognormal(fu, stwo ./ lambda);
  for (j in 1:n) {
    counts[, j] ~ neg_binomial_2_log(
      (log(size_factors[j]) + mu),
      1 ./ delta
    );
  }
}
