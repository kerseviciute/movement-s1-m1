import mne
import pandas as pd
import numpy as np

pd.to_pickle(snakemake, ".action_potential.py.pkl")
# snakemake = pd.read_pickle(".action_potential.py.pkl")

with open(f"{snakemake.scriptdir}/action_potential_methods.py", "r") as file:
    exec(file.read())

print(f"Detecting action potentials in cell {snakemake.wildcards['animal_id']} {snakemake.wildcards['cell_name']}")

vm = pd.read_pickle(snakemake.input["vm"])
data = vm.get_data()
time = vm.times

aps = []

for i, _ in enumerate(vm.ch_names):
    signal = data[i, :]
    differential = np.diff(signal)

    channel_ap = find_ap(differential, time, threshold = 1.5)

    if len(channel_ap) >= 1:
        channel_ap["Channel"] = i
        aps.append(channel_ap)

if len(aps) > 1:
    aps = pd.concat(aps)
else:
    aps = pd.DataFrame(
        columns = [ "EventStart", "EventEnd", "Start", "End", "Channel" ]
    )

# Save results
aps.to_csv(snakemake.output["action_potentials"], index = False)
