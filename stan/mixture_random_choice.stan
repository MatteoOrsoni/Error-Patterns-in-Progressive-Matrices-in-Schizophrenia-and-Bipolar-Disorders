// Mixture of a systematic process and random choice.
//
// Each response comes either from the systematic process (probability 1 - lambda)
// or from a random choice among the response options of the item (probability
// lambda), in which case each option is equally likely, so that each error type
// is chosen in proportion to its number of distractors on that item.
//
// include_correct = 0 (errors only, as Model 1): y in 1..K (R, D, IC, WP);
//   P(y = k) = (1 - lambda) * softmax(eta)[k] + lambda * q[i, k],
//   with q[i, k] = n_ik / sum_k n_ik (distractor composition of item i).
// include_correct = 1 (all answered items): y in 0..K, 0 = correct response;
//   the systematic process answers correctly with probability pc and otherwise
//   produces error type k with probability softmax(eta)[k]; random choice is
//   uniform over the 1 + D options (correct option + D distractors):
//   P(y = 0) = (1 - lambda) * pc + lambda / (1 + D),
//   P(y = k) = (1 - lambda) * (1 - pc) * softmax(eta)[k] + lambda * D / (1 + D) * q[i, k].
//
// Systematic error-type part: identical to multinomial_exchangeable.stan
//   eta[n, k] = alpha[k] + X[n] * beta[, k] + u[j, k] + v[i, k], R = reference.
// Random-choice share: logit(lambda[n]) = gamma[g] + w[j], with g the group of
//   participant j and w[j] ~ Normal(0, sigma_w) (participant variation).
// Accuracy part (include_correct = 1 only):
//   logit(pc[n]) = a0 + Xa[n] * b_acc + ua[j] + va[i].
// Salience (salient = k > 0): random choice favours the distractors of error
//   type k, each of which has weight exp(delta) instead of 1 (delta = 0: uniform
//   choice); salient = 0 gives the uniform random choice described above.
functions {
  // Probabilities of the random-choice component: element 1 = correct response,
  // elements 2..K+1 = error types.
  vector random_choice(row_vector qi, int D, int include_correct, int salient,
                       array[] real delta) {
    int K = cols(qi);
    if (salient == 0) {
      if (include_correct == 0) return append_row(0, qi');
      return append_row(1.0 / (1 + D), (1.0 * D / (1 + D)) * qi');
    }
    vector[K] wts = D * qi';
    wts[salient] *= exp(delta[1]);
    if (include_correct == 0) return append_row(0, wts / sum(wts));
    return append_row(1, wts) / (1 + sum(wts));
  }
}
data {
  int<lower=1> N;
  int<lower=2> K;                              // error types (R reference)
  int<lower=0> P;                              // predictors of error type (0 = none)
  int<lower=1> J;
  int<lower=1> I;
  int<lower=1> G;                              // groups with their own lambda
  array[N] int<lower=1, upper=J> subject;
  array[N] int<lower=1, upper=I> item;
  array[J] int<lower=1, upper=G> subject_group;
  int<lower=0, upper=1> include_correct;
  array[N] int<lower=0, upper=K> y;            // 0 = correct (include_correct = 1 only)
  matrix[N, P] X;
  int<lower=0> Pa;                             // predictors of accuracy
  matrix[N, Pa] Xa;
  matrix<lower=0, upper=1>[I, K] q;            // distractor composition, rows sum to 1
  int<lower=1> D;                              // distractors per item
  int<lower=0, upper=K> salient;               // error type favoured by random choice (0 = none)
}
transformed data {
  int K_minus_1 = K - 1;
  int A = include_correct;                     // size of the accuracy parameters
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

  vector[G] gamma;
  real<lower=0> sigma_w;
  vector[J] z_w;

  array[A] real a0;
  vector[A * Pa] b_acc;
  array[A] real<lower=0> sigma_ua;
  vector[A * J] z_ua;
  array[A] real<lower=0> sigma_va;
  vector[A * I] z_va;

  array[salient > 0] real delta;               // log weight of the salient distractor type
}
transformed parameters {
  matrix[P, K_minus_1] beta = rep_matrix(mu_beta, K_minus_1)
                              + rep_matrix(tau_beta, K_minus_1) .* z_beta;
  matrix[J, K_minus_1] u = z_u .* rep_matrix(sigma_u', J);
  matrix[I, K_minus_1] v = z_v .* rep_matrix(sigma_v', I);
  vector[J] w = sigma_w * z_w;
}
model {
  matrix[N, K_minus_1] eta_nr = rep_matrix(alpha', N) + X * beta + u[subject] + v[item];

  alpha    ~ normal(0, 1.5);
  mu_beta  ~ normal(0, 0.5);
  tau_beta ~ exponential(1);
  to_vector(z_beta) ~ std_normal();
  sigma_u  ~ student_t(3, 0, 1);
  sigma_v  ~ student_t(3, 0, 1);
  to_vector(z_u) ~ std_normal();
  to_vector(z_v) ~ std_normal();

  delta   ~ normal(0, 1);
  gamma   ~ logistic(0, 1);                    // uniform prior on lambda of a typical participant
  sigma_w ~ normal(0, 1);
  z_w     ~ std_normal();

  if (include_correct) {
    a0       ~ normal(0, 1.5);
    b_acc    ~ normal(0, 1);
    sigma_ua ~ student_t(3, 0, 1);
    sigma_va ~ student_t(3, 0, 1);
    z_ua     ~ std_normal();
    z_va     ~ std_normal();
  }

  for (n in 1:N) {
    real lambda_logit = gamma[subject_group[subject[n]]] + w[subject[n]];
    vector[K] p_sys = softmax(append_row(0, eta_nr[n]'));
    if (salient > 0) {
      vector[K + 1] p_rand = random_choice(q[item[n]], D, include_correct, salient, delta);
      real lp_sys;
      if (include_correct == 0) {
        lp_sys = log(p_sys[y[n]]);
      } else {
        real pc_logit = a0[1] + Xa[n] * b_acc
                        + sigma_ua[1] * z_ua[subject[n]] + sigma_va[1] * z_va[item[n]];
        lp_sys = y[n] == 0 ? log_inv_logit(pc_logit)
                           : log1m_inv_logit(pc_logit) + log(p_sys[y[n]]);
      }
      target += log_mix(inv_logit(lambda_logit), log(p_rand[y[n] + 1]), lp_sys);
    } else if (include_correct == 0) {
      target += log_mix(inv_logit(lambda_logit), log(q[item[n], y[n]]), log(p_sys[y[n]]));
    } else {
      real pc_logit = a0[1] + Xa[n] * b_acc
                      + sigma_ua[1] * z_ua[subject[n]] + sigma_va[1] * z_va[item[n]];
      if (y[n] == 0) {
        target += log_mix(inv_logit(lambda_logit), -log1p(D), log_inv_logit(pc_logit));
      } else {
        target += log_mix(inv_logit(lambda_logit),
                          log(D) - log1p(D) + log(q[item[n], y[n]]),
                          log1m_inv_logit(pc_logit) + log(p_sys[y[n]]));
      }
    }
  }
}
generated quantities {
  vector[N] log_lik;
  vector[G] lambda_typical = inv_logit(gamma);     // participant with w = 0
  vector[G] lambda_mean;                           // mean lambda of the participants of each group
  matrix[G, K + 1] expected_counts = rep_matrix(0, G, K + 1);  // col 1 = correct
  vector[G] random_errors = rep_vector(0, G);      // expected errors from the random component
  {
    matrix[N, K_minus_1] eta_nr = rep_matrix(alpha', N) + X * beta + u[subject] + v[item];
    vector[G] n_subj = rep_vector(0, G);
    lambda_mean = rep_vector(0, G);
    for (j in 1:J) {
      lambda_mean[subject_group[j]] += inv_logit(gamma[subject_group[j]] + w[j]);
      n_subj[subject_group[j]] += 1;
    }
    lambda_mean = lambda_mean ./ n_subj;

    for (n in 1:N) {
      int g = subject_group[subject[n]];
      real lam = inv_logit(gamma[g] + w[subject[n]]);
      vector[K] p_sys = softmax(append_row(0, eta_nr[n]'));
      vector[K + 1] p_rand = random_choice(q[item[n]], D, include_correct, salient, delta);
      vector[K + 1] p_syst;
      if (include_correct == 0) {
        p_syst = append_row(0, p_sys);
      } else {
        real pc = inv_logit(a0[1] + Xa[n] * b_acc
                            + sigma_ua[1] * z_ua[subject[n]] + sigma_va[1] * z_va[item[n]]);
        p_syst = append_row(pc, (1 - pc) * p_sys);
      }
      vector[K + 1] p = lam * p_rand + (1 - lam) * p_syst;
      log_lik[n] = log(p[y[n] + 1]);
      expected_counts[g] += p';
      if (y[n] > 0) random_errors[g] += lam * p_rand[y[n] + 1] / p[y[n] + 1];
    }
  }
}
