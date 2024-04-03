import mne
import numpy as np
import pandas as pd

pd.to_pickle(snakemake, ".filter.py.pkl")
# snakemake = pd.read_pickle(".filter.py.pkl")

with open(f"{snakemake.scriptdir}/to_mne.py", "r") as file:
    exec(file.read())

sfreq = snakemake.params["sampling_rate"]
ch_type = snakemake.params["ch_type"]

raw = read_csv_as_mne(snakemake.input["raw"], sfreq = sfreq, ch_type = ch_type)

l_freq = snakemake.params["drop_below"]
h_freq = snakemake.params["drop_above"]

raw.filter(l_freq = l_freq, h_freq = h_freq, fir_design = "firwin", picks = ch_type)

if snakemake.params["scale"]:
    print("Scaling the channel data")
    data = raw.get_data().transpose()
    data = (data - np.mean(data)) / np.std(data)
    raw._data = data.transpose()
else:
    print("Channels were not scaled")

data = pd.DataFrame(raw.get_data())
data.index = raw.ch_names

data.to_csv(snakemake.output["filter"])
