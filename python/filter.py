import mne
import numpy as np
import pandas as pd

pd.to_pickle(snakemake, ".filter.py.pkl")
# snakemake = pd.read_pickle(".filter.py.pkl")

raw = pd.read_pickle(snakemake.input["raw"])

l_freq = snakemake.params["drop_below"]
h_freq = snakemake.params["drop_above"]

ch_type = snakemake.params["ch_type"]
raw.filter(l_freq = l_freq, h_freq = h_freq, fir_design = "firwin", picks = ch_type)

if snakemake.params["scale"]:
    print("Scaling the channel data")
    data = raw.get_data().transpose()
    data = (data - np.mean(data)) / np.std(data)
    raw._data = data.transpose()
else:
    print("Channels were not scaled")

pd.to_pickle(raw, snakemake.output["filter"])
