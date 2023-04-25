functions {
#include ode_pratt_fix.stan
}
transformed data {
  int<lower=1> N_timepoint = 1000;
  vector[3] breaks = [0, 300, 600]';
  real k_tg = 0.00185859360751066;
  real k_glucose = 0.139081199343437;
  real d_tg = 180;
  real d_glucose = 45;
  real d_glucose_bad = 180;
}
generated quantities {
  vector[N_timepoint] diet_tg;
  vector[N_timepoint] diet_glucose;
  vector[N_timepoint] diet_glucose_bad;
  for (t in 1:N_timepoint){
    diet_tg[t] = get_diet(t, breaks, k_tg, d_tg);
    diet_glucose[t] = get_diet(t, breaks, k_glucose, d_glucose);
    diet_glucose_bad[t] = get_diet(t, breaks, k_glucose, d_glucose_bad);
  }
}
