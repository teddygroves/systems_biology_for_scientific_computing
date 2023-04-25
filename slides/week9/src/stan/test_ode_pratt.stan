functions{
#include ode_pratt.stan
}
data {
  int<lower=1> N_species;
  int<lower=1> N_parameter;
  int<lower=1> N_reaction;
  int<lower=1> N_compartment;
  int<lower=1> N_timepoint;
  matrix[N_species,N_reaction] S_rxn;
  matrix[N_species,N_compartment] S_cpt;
  vector[N_compartment] compartment_volume;
  vector[N_species] initial_concentration;
  array[N_timepoint] real timepoint;
  real rel_tol;
  real abs_tol;
  int<lower=0> max_num_steps;
  vector<lower=0>[N_parameter] p;
}
generated quantities {
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
