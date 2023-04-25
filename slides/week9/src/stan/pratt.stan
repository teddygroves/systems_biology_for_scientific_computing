functions {
#include ode_pratt.stan
}
data {
  int<lower=1> N_species;
  int<lower=1> N_reaction;
  int<lower=1> N_compartment;
  int<lower=1> N_parameter;
  int<lower=1> N_fixed_parameter;
  int<lower=1> N_timepoint;
  int<lower=0,upper=N_timepoint*N_species> N_measurement;
  // network information
  matrix[N_species,N_reaction] S_rxn;
  matrix[N_species,N_compartment] S_cpt;
  // known quantities
  vector[N_compartment] compartment_volume;
  vector[N_species] initial_concentration;
  array[N_timepoint] real timepoint;
  vector[N_fixed_parameter] fixed_parameter_value;
  array[N_fixed_parameter] int<lower=1,upper=N_parameter> fixed_parameter_ix;
  array[N_parameter-N_fixed_parameter] int<lower=1,upper=N_parameter> free_parameter_ix;
  // measurement information
  array[N_measurement] int<lower=1,upper=N_timepoint> y_timepoint_ix;
  array[N_measurement] int<lower=1,upper=N_species> y_species_ix;
  vector[N_measurement] y;
  vector[N_measurement] y_sd;
  // priors
  array[2] vector[N_parameter-N_fixed_parameter] prior_p;
  // run configuration
  int<lower=0,upper=1> likelihood;
  real rel_tol;
  real abs_tol;
  int<lower=0> max_num_steps;
}
parameters {
  vector<lower=0>[N_parameter-N_fixed_parameter] p_free;
}
transformed parameters {
  vector[N_parameter] p;
  p[free_parameter_ix] = p_free;
  p[fixed_parameter_ix] = fixed_parameter_value;
  array[N_timepoint] vector[N_species] conc = 
    ode_bdf_tol(get_dxdt_pratt,
                initial_concentration, 
                0, 
                timepoint, 
                rel_tol,
                abs_tol,
                max_num_steps,
                p, 
                S_rxn, 
                S_cpt, 
                compartment_volume);
  array[N_timepoint] vector[N_reaction] flux;
  for (t in 1:N_timepoint){
    flux[t] = get_flux_pratt(t, conc[t], p);
  }
}
model {
  // prior model
  p_free ~ lognormal(prior_p[1], prior_p[2]);
  // measurement model
  if (likelihood){
    for (m in 1:N_measurement){
      y[m] ~ lognormal(log(conc[y_timepoint_ix[m]][y_species_ix[m]]), y_sd[m]);
    }
  }
}
generated quantities {
  vector[N_measurement] yrep;
  for (m in 1:N_measurement){
    yrep[m] = lognormal_rng(log(conc[y_timepoint_ix[m]][y_species_ix[m]]), y_sd[m]);
  }
} 
