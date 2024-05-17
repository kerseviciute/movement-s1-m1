import mne
import pandas as pd
import numpy as np

pd.to_pickle(snakemake, ".movement_offset.py.pkl")
# snakemake = pd.read_pickle(".movement_offset.py.pkl")

with open(f"{snakemake.scriptdir}/to_mne.py", "r") as file:
    exec(file.read())

movement = pd.read_csv(snakemake.input["movement"])
sfreq = snakemake.params["sfreq"]
raw = read_csv_as_mne(snakemake.input["raw"], sfreq = sfreq)

# Evaluate channel type based on input data
if snakemake.wildcards["type"] == "emg":
    # High-pass filter the data at 2 Hz
    raw.filter(l_freq = 2, h_freq = None, fir_design = "firwin", picks = "emg")

sample_data = raw.get_data()

episode_data = []

for i, episode in movement.iterrows():
    start = episode["EventEnd"] - round(sfreq * 0.5)
    end = episode["EventEnd"] + round(sfreq * 0.5)
    channel = episode["Channel"]
    episode_data.append(sample_data[channel][start:end])

episode_data = np.array(episode_data)
episode_data = pd.DataFrame(episode_data)

episode_data.index = movement["ID"]

episode_data.to_csv(snakemake.output["offset"])
