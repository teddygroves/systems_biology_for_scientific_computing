functions {
  vector dxdt(real t, vector x, real ka, real vm, real v, real km){
    real c = x[2] / v;
    return [
      -ka * x[1],
      ka * x[1] - vm * c / (km + c)
    ]';
  }
}
data {
  int<lower=1> N_timepoint;
  int<lower=0,upper=N_timepoint> N_measurement;
  array[N_timepoint] real timepoint;
  vector[2] initial_concentration;
  // measurement information
  array[N_measurement] int<lower=1,upper=N_timepoint> y_timepoint_ix;
  array[N_measurement] vector[2] y;
  // run configuration
  int<lower=0,upper=1> likelihood;
}
parameters {
  real<lower=0> sigma;
  real<lower=0> ka;
  real<lower=0> vm;
  real<lower=0> v;
  real<lower=0> km;
}
transformed parameters {
array[N_timepoint] vector[2] conc = ode_bdf(dxdt,
                                            initial_concentration, 
                                            0, 
                                            timepoint, 
                                            ka,
                                            vm,
                                            v,
                                            km); 
}
model {
  sigma ~ normal(0, 1);
  ka ~ lognormal(log(2.5), 3);
  vm ~ lognormal(log(2.5), 3);
  v ~ lognormal(log(35), 0.5);
  km ~ lognormal(log(10), 0.5);
  if (likelihood){
    for (m in 1:N_measurement){
      for (s in 1:2){
        y[m][s] ~ lognormal(log(conc[y_timepoint_ix[m]][s]), sigma);
      }
    }
  }
}
generated quantities {
  array[N_measurement] vector[2] yrep;
    for (m in 1:N_measurement){
      for (s in 1:2)
        yrep[m,s] = lognormal_rng(log(conc[y_timepoint_ix[m]][s]), sigma);
    }
}
