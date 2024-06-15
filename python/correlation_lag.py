import mne
import pandas as pd
import numpy as np

pd.to_pickle(snakemake, ".correlation_lag.py.pkl")
# snakemake = pd.read_pickle(".correlation_lag.py.pkl")

with open(f"{snakemake.scriptdir}/to_mne.py", "r") as file:
    exec(file.read())

with open(f"{snakemake.scriptdir}/correlation_lag_methods.py", "r") as file:
    exec(file.read())

sfreq = snakemake.params["sfreq"]
emg = read_csv_as_mne(snakemake.input["emg"], sfreq = sfreq)
vm = read_csv_as_mne(snakemake.input["vm"], sfreq = sfreq)

print(f"Calculating correlation in cell {snakemake.wildcards['animal_id']} {snakemake.wildcards['cell_name']}")

lags, correlations = correlation_lag(emg, vm, lag_range = [-20_000, 20_000])

# Save results
correlations.to_csv(snakemake.output["correlation"], header = False, index = False)
