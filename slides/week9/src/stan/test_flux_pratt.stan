functions {
#include ode_pratt.stan
}
data {
  int<lower=1> N_species;
  int<lower=1> N_reaction;
  int<lower=1> N_parameter;
  vector<lower=0>[N_species] x;
  vector<lower=0>[N_parameter] p;
  real<lower=0> t;
}
generated quantities {
  vector[N_reaction] flux = get_flux_pratt(t, x, p);
}
