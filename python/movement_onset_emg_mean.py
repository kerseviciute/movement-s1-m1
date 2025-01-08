import mne
import pandas as pd
import numpy as np

pd.to_pickle(snakemake, ".movement_onset_emg_mean.py.pkl")
# snakemake = pd.read_pickle(".movement_onset_emg_mean.py.pkl")

with open(f"{snakemake.scriptdir}/to_mne.py", "r") as file:
    exec(file.read())

with open(f"{snakemake.scriptdir}/correlation_lag_methods.py", "r") as file:
    exec(file.read())

sfreq = snakemake.params["sfreq"]
emg = read_csv_as_mne(snakemake.input["emg"], sfreq = sfreq)

onset_time = (np.max(emg.times) + 1 / sfreq) / 2
time_window = 0.4

res = []
for i, episode_id in enumerate(emg.ch_names):
    emg_data = emg.copy().crop(tmin = onset_time - time_window, tmax = onset_time + time_window).pick(picks = [i]).get_data()[ 0, : ]
    filtered = filter_emg(emg_data)

    onset_idx = round(len(emg_data) / 2)

    before = np.mean(filtered[ :onset_idx ])
    after = np.mean(filtered[ onset_idx: ])
    diff = after - before

    res.append(pd.DataFrame({
        "ID": [episode_id],
        "EMG": [diff]
    }))

res = pd.concat(res)
res.reset_index(drop = True, inplace = True)

# Save results
res.to_csv(snakemake.output["mean"], index = False)
