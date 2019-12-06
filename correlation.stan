data {
  int<lower=1> N;
  vector[2] x[N];
}
parameters {
  vector[2] mu;
  real<lower=0> sigma[2];
  real<lower=-1, upper=1> rho;
}
transformed parameters {
  cov_matrix[2] cov;
  cov = [[sigma[1]^2,               sigma[1] * sigma[2] * rho],
        [sigma[1] * sigma[2] * rho, sigma[2]^2              ]];
}
model {
  x ~ multi_normal(mu, cov);
  sigma ~ student_t(3, 0, 2);
  mu ~ normal(0, 1);
}
