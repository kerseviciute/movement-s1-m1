import mne
import numpy as np
import pandas as pd

pd.to_pickle(snakemake, ".vm_statistics.py.pkl")
# snakemake = pd.read_pickle(".vm_statistics.py.pkl")


vm_data = pd.read_csv(snakemake.input["vm"], index_col = 0).to_numpy()

action_potentials = pd.read_csv(snakemake.input["action_potentials"])
movement = pd.read_csv(snakemake.input["movement"])
rest = pd.read_csv(snakemake.input["rest"])

# Movement

stats = []

for i, episode in movement.iterrows():
    start = episode["EventStart"]
    end = episode["EventEnd"]
    channel = episode["Channel"]
    eid = episode["ID"]

    signal_range = np.arange(start, end)

    channel_aps = action_potentials[action_potentials["Channel"] == channel]
    episode_aps = np.logical_and(
        channel_aps["EventStart"] >= start,
        channel_aps["EventStart"] <= end
    )

    n_ap = sum(episode_aps)

    if sum(episode_aps) >= 1:
        # Get APs within the episode and remove their ranges
        episode_aps = channel_aps[episode_aps]
        for j, ap in episode_aps.iterrows():
            ap_start = round(ap["EventStart"])
            ap_end = round(ap["EventEnd"]) + 50

            ap_range = np.arange(ap_start, ap_end)
            signal_range = [i for i in signal_range if i not in ap_range]

    signal = vm_data[channel, start:end]
    signal_no_ap = vm_data[channel, signal_range]

    stats.append(pd.DataFrame({
        "ID": [eid],
        "Mean": [np.mean(signal)],
        "SD": [np.std(signal)],
        "MeanNoAP": [np.mean(signal_no_ap)],
        "SDNoAP": [np.std(signal_no_ap)],
        "NumberOfAP": [n_ap]
    }))

# Rest

for i, episode in rest.iterrows():
    start = episode["EventStart"]
    end = episode["EventEnd"]
    channel = episode["Channel"]
    eid = episode["ID"]

    signal_range = np.arange(start, end)

    channel_aps = action_potentials[action_potentials["Channel"] == channel]
    episode_aps = np.logical_and(
        channel_aps["EventStart"] >= start,
        channel_aps["EventStart"] <= end
    )

    n_ap = sum(episode_aps)

    if sum(episode_aps) >= 1:
        # Get APs within the episode and remove their ranges
        episode_aps = channel_aps[episode_aps]
        for j, ap in episode_aps.iterrows():
            ap_start = round(ap["EventStart"])
            ap_end = round(ap["EventEnd"]) + 50

            ap_range = np.arange(ap_start, ap_end)
            signal_range = [i for i in signal_range if i not in ap_range]

    signal = vm_data[channel, start:end]
    signal_no_ap = vm_data[channel, signal_range]

    stats.append(pd.DataFrame({
        "ID": [eid],
        "Mean": [np.mean(signal)],
        "SD": [np.std(signal)],
        "MeanNoAP": [np.mean(signal_no_ap)],
        "SDNoAP": [np.std(signal_no_ap)],
        "NumberOfAP": [n_ap]
    }))

stats = pd.concat(stats)
stats = stats.reset_index(drop = True)

stats.to_csv(snakemake.output["statistics"], index = False)
