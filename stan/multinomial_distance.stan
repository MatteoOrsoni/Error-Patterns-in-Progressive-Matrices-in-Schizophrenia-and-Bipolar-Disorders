// Model 1 at the level of the single distractor, with distractor distance.
//
// Outcome: the distractor chosen on an item, given that an error occurred.
// Each item i offers n_opt[i] distractors; distractor d has an error category
// opt_cat[i, d] (1 = R, reference; 2..K = D, IC, WP) and a distance from the
// correct response opt_dist[i, d] (number of departures, consensus coding).
// Utility of distractor d for response n (participant j, item i):
//   V[n, d] = eta[n, opt_cat[i, d]] + (gamma + X[n] * delta) * opt_dist[i, d]
//   eta[n, 1] = 0;  eta[n, k] = alpha[k] + X[n] * beta[, k] + u[j, k] + v[i, k]  (k > 1)
// P(d) = softmax(V[n, 1:n_opt[i]]). A response label can identify more than
// one distractor (the "r.ic" label in three items, see the R script): the
// likelihood of a response is the sum of the probabilities of the distractors
// flagged in chosen[n, ].
// With use_distance = 0 and use_interaction = 0 the model is Model 1 with the
// number of distractors of each category on each item as an offset: the
// probability of an error category is n_ik * exp(eta_k) / sum_k' n_ik' * exp(eta_k').
// The group coefficients beta have the same meaning as in Model 1; with distance,
// they compare a distractor of category k with an R distractor at the same
// distance from the correct response.
// Priors and the structure of alpha, beta, u and v are those of
// multinomial_exchangeable.stan. log_lik is the log-likelihood of the chosen
// distractor, for comparisons among these models; log_lik_cat that of the error
// category (the outcome of Model 1), for the comparison with Model 1.
data {
  int<lower=1> N;
  int<lower=2> K;
  int<lower=1> P;
  int<lower=1> J;
  int<lower=1> I;
  int<lower=1> D_max;
  array[N] int<lower=1, upper=J> subject;
  array[N] int<lower=1, upper=I> item;
  matrix[N, P] X;
  array[I] int<lower=1, upper=D_max> n_opt;
  array[I, D_max] int<lower=0, upper=K> opt_cat;     // 0 = padding
  array[I, D_max] real opt_dist;
  array[N, D_max] int<lower=0, upper=1> chosen;
  int<lower=0, upper=1> use_distance;                // main effect of distance
  int<lower=0, upper=1> use_interaction;             // predictors x distance
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
  array[use_distance] real gamma;                    // log-odds per unit of distance (reference group)
  vector[use_interaction ? P : 0] delta;             // change in the distance slope per predictor
}
transformed parameters {
  matrix[P, K_minus_1] beta = rep_matrix(mu_beta, K_minus_1)
                              + rep_matrix(tau_beta, K_minus_1) .* z_beta;
  matrix[J, K_minus_1] u = z_u .* rep_matrix(sigma_u', J);
  matrix[I, K_minus_1] v = z_v .* rep_matrix(sigma_v', I);
}
model {
  alpha    ~ normal(0, 1.5);
  mu_beta  ~ normal(0, 0.5);
  tau_beta ~ exponential(1);
  to_vector(z_beta) ~ std_normal();
  sigma_u  ~ student_t(3, 0, 1);
  sigma_v  ~ student_t(3, 0, 1);
  to_vector(z_u) ~ std_normal();
  to_vector(z_v) ~ std_normal();
  gamma    ~ normal(0, 1);
  delta    ~ normal(0, 0.5);
  {
    matrix[N, K_minus_1] eta_nr = rep_matrix(alpha', N) + X * beta + u[subject] + v[item];
    for (n in 1:N) {
      int i = item[n];
      int m = n_opt[i];
      vector[K] eta = append_row(0, eta_nr[n]');
      real slope = 0;
      vector[m] V;
      if (use_distance) slope += gamma[1];
      if (use_interaction) slope += X[n] * delta;
      for (d in 1:m) V[d] = eta[opt_cat[i, d]] + slope * opt_dist[i, d];
      {
        real lse_all = log_sum_exp(V);
        vector[m] V_chosen = V;
        for (d in 1:m) if (chosen[n, d] == 0) V_chosen[d] = negative_infinity();
        target += log_sum_exp(V_chosen) - lse_all;
      }
    }
  }
}
generated quantities {
  vector[N] log_lik;          // pointwise log-likelihood of the chosen distractor (PSIS-LOO)
  vector[N] log_lik_cat;      // of the error category (sum over the distractors of the chosen
                              // category): same outcome as Model 1, for the comparison with it
  {
    matrix[N, K_minus_1] eta_nr = rep_matrix(alpha', N) + X * beta + u[subject] + v[item];
    for (n in 1:N) {
      int i = item[n];
      int m = n_opt[i];
      vector[K] eta = append_row(0, eta_nr[n]');
      real slope = 0;
      vector[m] V;
      if (use_distance) slope += gamma[1];
      if (use_interaction) slope += X[n] * delta;
      for (d in 1:m) V[d] = eta[opt_cat[i, d]] + slope * opt_dist[i, d];
      {
        vector[m] V_chosen = V;
        for (d in 1:m) if (chosen[n, d] == 0) V_chosen[d] = negative_infinity();
        log_lik[n] = log_sum_exp(V_chosen) - log_sum_exp(V);
      }
      {
        int c = 0;
        vector[m] V_cat = V;
        for (d in 1:m) if (chosen[n, d] == 1) c = opt_cat[i, d];
        for (d in 1:m) if (opt_cat[i, d] != c) V_cat[d] = negative_infinity();
        log_lik_cat[n] = log_sum_exp(V_cat) - log_sum_exp(V);
      }
    }
  }
}
