functions {
  real mm_one_substrate(real vmax, real km, real sub){
    return vmax * sub / (km + sub);
  }
  vector get_flux_simple(real t, vector x, vector p){
    return [mm_one_substrate(p[1], p[2], x[1])]';
  }
  vector get_dxdt_simple(real t, 
                         vector x, 
                         vector p, 
                         matrix S_rxn,
                         matrix S_cpt,
                         vector v){
    return (S_rxn * get_flux_simple(t, p, x)) ./ (S_cpt * v);
  }
}
data {
  // dimensions
  int<lower=1> N_species;
  int<lower=1> N_reaction;
  int<lower=1> N_compartment;
  int<lower=1> N_timepoint;
  int<lower=0,upper=N_timepoint> N_measurement;
  // network information
  matrix[N_species,N_reaction] S_rxn;
  matrix[N_species,N_compartment] S_cpt;
  // known quantities
  vector[N_compartment] compartment_volume;
  vector[N_species] initial_concentration;
  array[N_timepoint] real timepoint;
  // measurement information
  array[N_measurement] int<lower=1,upper=N_timepoint> y_timepoint_ix;
  array[N_measurement] int<lower=1,upper=N_species> y_species_ix;
  vector[N_measurement] y;
  vector[N_measurement] y_sd;
  // priors
  array[2] vector[2] prior_p;
  // run configuration
  int<lower=0,upper=1> likelihood;
}
parameters {
  vector<lower=0>[2] p;
}
transformed parameters {
  // species concentrations at measurement times
  array[N_timepoint] vector[N_species] conc = 
    ode_bdf(get_dxdt_simple,
            initial_concentration, 
            0, 
            timepoint, 
            p, 
            S_rxn, 
            S_cpt, 
            compartment_volume);
}
model {
  // prior model
  p ~ lognormal(prior_p[1], prior_p[2]);
  // measurement model
  if (likelihood){
    for (m in 1:N_measurement){
      y[m] ~ lognormal(log(conc[y_timepoint_ix[m]][y_species_ix[m]]), y_sd[m]);
    }
  }
}
generated quantities {
  // simulated measurements
  vector[N_measurement] yrep;
  for (m in 1:N_measurement){
    yrep[m] = lognormal_rng(log(conc[y_timepoint_ix[m]][y_species_ix[m]]), 
                                y_sd[m]);
  }
} 
