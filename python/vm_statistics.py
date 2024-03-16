import mne
import numpy as np
import pandas as pd

pd.to_pickle(snakemake, ".vm_statistics.py.pkl")
# snakemake = pd.read_pickle(".vm_statistics.py.pkl")


vm = pd.read_pickle(snakemake.input["vm"])
vm_data = vm.get_data()

action_potentials = pd.read_csv(snakemake.input["action_potentials"])
movement = pd.read_csv(snakemake.input["movement"])
rest = pd.read_csv(snakemake.input["rest"])

#
# Evaluate including action potentials
#

# Movement

movement_vm = []

for i, episode in movement.iterrows():
    start = episode["EventStart"]
    end = episode["EventEnd"]
    channel = episode["Channel"]

    signal = vm_data[channel, start:end]
    movement_vm.append(signal)

movement_vm = np.concatenate(movement_vm)

# Rest

rest_vm = []

for i, episode in rest.iterrows():
    start = episode["EventStart"]
    end = episode["EventEnd"]
    channel = episode["Channel"]

    signal = vm_data[channel, start:end]
    rest_vm.append(signal)

rest_vm = np.concatenate(rest_vm)

#
# Evaluate excluding action potentials
#

# Movement

movement_vm_no_ap = []
movement_n_ap = 0

for i, episode in movement.iterrows():
    start = episode["EventStart"]
    end = episode["EventEnd"]
    channel = episode["Channel"]

    signal_range = np.arange(start, end)

    channel_aps = action_potentials[action_potentials["Channel"] == channel]
    episode_aps = np.logical_and(
        channel_aps["EventStart"] >= start,
        channel_aps["EventStart"] <= end
    )

    movement_n_ap += sum(episode_aps)

    if sum(episode_aps) >= 1:
        # Get APs within the episode and remove their ranges
        episode_aps = channel_aps[episode_aps]
        for j, ap in episode_aps.iterrows():
            ap_start = int(ap["EventStart"])
            ap_end = int(ap["EventEnd"])

            ap_range = np.arange(ap_start, ap_end)
            signal_range = [i for i in signal_range if i not in ap_range]

    signal = vm_data[channel, signal_range]
    movement_vm_no_ap.append(signal)

movement_vm_no_ap = np.concatenate(movement_vm_no_ap)


# Rest

rest_vm_no_ap = []
rest_n_ap = 0

for i, episode in rest.iterrows():
    start = episode["EventStart"]
    end = episode["EventEnd"]
    channel = episode["Channel"]

    signal_range = np.arange(start, end)

    channel_aps = action_potentials[action_potentials["Channel"] == channel]
    episode_aps = np.logical_and(
        channel_aps["EventStart"] >= start,
        channel_aps["EventStart"] <= end
    )

    rest_n_ap += sum(episode_aps)

    if sum(episode_aps) >= 1:
        # Get APs within the episode and remove their ranges
        episode_aps = channel_aps[episode_aps]
        for j, ap in episode_aps.iterrows():
            ap_start = int(ap["EventStart"])
            ap_end = int(ap["EventEnd"])

            ap_range = np.arange(ap_start, ap_end)
            signal_range = [i for i in signal_range if i not in ap_range]

    signal = vm_data[channel, signal_range]
    rest_vm_no_ap.append(signal)

rest_vm_no_ap = np.concatenate(rest_vm_no_ap)


statistics = pd.DataFrame({
    "Mean": [np.mean(movement_vm), np.mean(rest_vm)],
    "MeanNoAP": [np.mean(movement_vm_no_ap), np.mean(rest_vm_no_ap)],
    "SD": [np.std(movement_vm), np.std(rest_vm)],
    "SDNoAP": [np.std(movement_vm_no_ap), np.std(rest_vm_no_ap)],
    "NumberOfAP": [movement_n_ap, rest_n_ap],
    "Episode": ["Movement", "Rest"]
})

statistics.to_csv(snakemake.output["statistics"], index = False)
