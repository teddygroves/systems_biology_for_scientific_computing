import cmdstanpy
import numpy as np
import pandas as pd

N_TIMEPOINT = 50
N_MEASUREMENT = 10
DURATION = 1
INITIAL_CONCENTRATION = [10, 0.1]
TRUTH = {
    "sigma": 0.1,
    "ka": 3,
    "vm": 2,
    "v": 50,
    "km": 18,
}

def main():
    model = cmdstanpy.CmdStanModel(stan_file="model.stan")
    timepoint = np.linspace(0, DURATION, N_TIMEPOINT + 1)[1:]  # miss out t0
    measurement_ix = np.sort(
        np.random.choice( N_TIMEPOINT, size=N_MEASUREMENT, replace=False)
    )
    y = np.ones((N_MEASUREMENT, 2))
    data = {
        "N_timepoint": len(timepoint),
        "N_measurement": len(measurement_ix),
        "timepoint": timepoint,
        "initial_concentration": INITIAL_CONCENTRATION,
        "y_timepoint_ix": measurement_ix + 1,
        "y": y,
        "likelihood": 0
    }
    mcmc = model.sample(
        data=data,
        chains=1,
        fixed_param=True,
        inits=TRUTH,
        iter_sampling=1
    )
    sim = pd.DataFrame(mcmc.stan_variable("yrep")[0], columns=["a", "b"])
    sim.index.name = "measurement"
    sim["time"] = timepoint[measurement_ix]
    sim.to_csv("simulated_measurements.csv")

if __name__ == "__main__":
    main()


