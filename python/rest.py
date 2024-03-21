import mne
import numpy as np
import pandas as pd

pd.to_pickle(snakemake, ".rest.py.pkl")
# snakemake = pd.read_pickle(".rest.py.pkl")

with open(f"{snakemake.scriptdir}/movement_methods.py", "r") as file:
    exec(file.read())

with open(f"{snakemake.scriptdir}/rest_methods.py", "r") as file:
    exec(file.read())

data = pd.read_pickle(snakemake.input["data"])
sample_data = data.get_data()

min_break = data.info["sfreq"] * snakemake.params["maxTimeApart"]
min_event_length = data.info["sfreq"] * snakemake.params["minEventLength"]
percentile = snakemake.params["percentile"]

# Apply low frequency filter on all channels
filtered_data = []
for channel in sample_data:
    filtered_data.append(filter_emg(channel, low_freq = 20, sfreq = data.info["sfreq"]))
filtered_data = np.array(filtered_data)

threshold = np.percentile(np.abs(sample_data), percentile)

episodes = []

for i, _ in enumerate(data.ch_names):
    signal = filtered_data[i]

    ch_episodes = detect_rest_episodes(
        signal = signal,
        threshold = threshold,
        sfreq = data.info["sfreq"],
        min_break = min_break
    )
    ch_episodes["Channel"] = i
    episodes.append(ch_episodes)

episodes = pd.concat(episodes)
episodes = episodes[ episodes["EventLength"] > min_event_length ]

episodes = episodes.reset_index(drop = True)

print(f"Number of detected rest episodes: {len(episodes)}")

episodes.to_csv(snakemake.output["episodes"], index = False)
