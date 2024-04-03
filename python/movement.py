import numpy as np
import pandas as pd

pd.to_pickle(snakemake, ".movement.py.pkl")
# snakemake = pd.read_pickle(".movement.py.pkl")

with open(f"{snakemake.scriptdir}/movement_methods.py", "r") as file:
    exec(file.read())

sample_data = pd.read_csv(snakemake.input["data"], index_col = 0)
trials = sample_data.index

sfreq = snakemake.params["sfreq"]
min_break = sfreq * snakemake.params["maxTimeApart"]
min_event_length = sfreq * snakemake.params["minEventLength"]
percentile = snakemake.params["percentile"]

filtered_data = []
for i, _ in enumerate(trials):
    trial = sample_data.iloc[i]
    filtered_data.append(filter_emg(trial, low_freq = 20, sfreq = sfreq))
filtered_data = np.array(filtered_data)

threshold = np.percentile(np.abs(sample_data), percentile)

episodes = []

for i, _ in enumerate(trials):
    signal = filtered_data[i]

    ch_episodes = detect_movement_episodes(
        signal = signal,
        threshold = threshold,
        sfreq = sfreq,
        min_break = min_break
    )
    ch_episodes["Channel"] = i
    episodes.append(ch_episodes)

episodes = pd.concat(episodes)
episodes = episodes[ episodes["EventLength"] > min_event_length ]

episodes = episodes.reset_index(drop = True)
episodes["ID"] = [ "M" + str(index) for index in episodes.index ]

print(f"Number of detected movement episodes: {len(episodes)}")

episodes.to_csv(snakemake.output["episodes"], index = False)
