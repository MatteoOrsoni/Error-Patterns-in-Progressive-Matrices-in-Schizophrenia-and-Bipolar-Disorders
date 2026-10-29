// Bayesian multinomial (categorical-logit) mixed-effects model for error types.
//
// Outcome: type of error on an item, given that an error occurred.
//   Category 1 = repetition (R, reference); categories 2..K = D, IC, WP.
// Linear predictor for response n (participant j, item i) and category k > 1:
//   eta[n, k] = alpha[k] + X[n] * beta[, k] + u[j, k] + v[i, k]
// - beta[p, k]: coefficients of predictor p, partially pooled across the K-1
//   contrasts (exchangeable prior): beta[p, k] = mu_beta[p] + tau_beta[p] * z_beta[p, k]
// - u, v: crossed random intercepts for participants and items (non-centred),
//   one set per contrast, independent across contrasts.
data {
  int<lower=1> N;                              // number of error responses
  int<lower=2> K;                              // number of categories (incl. reference)
  int<lower=1> P;                              // number of predictors
  int<lower=1> J;                              // number of participants
  int<lower=1> I;                              // number of items
  array[N] int<lower=1, upper=J> subject;
  array[N] int<lower=1, upper=I> item;
  array[N] int<lower=1, upper=K> y;
  matrix[N, P] X;
}
transformed data {
  int K_minus_1 = K - 1;
}
parameters {
  vector[K_minus_1] alpha;
  vector[P] mu_beta;
  vector<lower=0>[P] tau_beta;
  matrix[P, K_minus_1] z_beta;
  vector<lower=0>[K_minus_1] sigma_u;
  matrix[J, K_minus_1] z_u;
  vector<lower=0>[K_minus_1] sigma_v;
  matrix[I, K_minus_1] z_v;
}
transformed parameters {
  matrix[P, K_minus_1] beta = rep_matrix(mu_beta, K_minus_1)
                              + rep_matrix(tau_beta, K_minus_1) .* z_beta;
  matrix[J, K_minus_1] u = z_u .* rep_matrix(sigma_u', J);
  matrix[I, K_minus_1] v = z_v .* rep_matrix(sigma_v', I);
}
model {
  matrix[N, K_minus_1] eta_nr = rep_matrix(alpha', N) + X * beta + u[subject] + v[item];

  alpha    ~ normal(0, 1.5);
  mu_beta  ~ normal(0, 0.5);
  tau_beta ~ exponential(1);
  to_vector(z_beta) ~ std_normal();
  sigma_u  ~ student_t(3, 0, 1);    // half-t because of the lower bound
  sigma_v  ~ student_t(3, 0, 1);
  to_vector(z_u) ~ std_normal();
  to_vector(z_v) ~ std_normal();

  for (n in 1:N)
    y[n] ~ categorical_logit(append_row(0, eta_nr[n]'));
}
generated quantities {
  // Pointwise log-likelihood and posterior predictive draws. Not used in the
  // paper's analyses; kept unchanged because the program must be identical to
  // the one used for the paper to reproduce its posterior draws exactly
  // (removing the random-number calls could alter the sampler's random stream).
  vector[N] log_lik;
  array[N] int y_rep;
  {
    matrix[N, K_minus_1] eta_nr = rep_matrix(alpha', N) + X * beta + u[subject] + v[item];
    for (n in 1:N) {
      vector[K] eta = append_row(0, eta_nr[n]');
      log_lik[n] = categorical_logit_lpmf(y[n] | eta);
      y_rep[n]   = categorical_logit_rng(eta);
    }
  }
}
