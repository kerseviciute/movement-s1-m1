import mne
import pandas as pd
import numpy as np

pd.to_pickle(snakemake, ".correlation_lag.py.pkl")
# snakemake = pd.read_pickle(".correlation_lag.py.pkl")

with open(f"{snakemake.scriptdir}/correlation_lag_methods.py", "r") as file:
    exec(file.read())

emg = pd.read_pickle(snakemake.input["emg"])
vm = pd.read_pickle(snakemake.input["vm"])

print(f"Calculating correlation in cell {snakemake.wildcards['animal_id']} {snakemake.wildcards['cell_name']}")

lags, correlations = correlation_lag(emg, vm, lag_range = [-20_000, 20_000])

# Save results
correlations.to_csv(snakemake.output["correlation"], header = False, index = False)
