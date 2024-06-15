import numpy as np
import pandas as pd
from numpy.fft import rfft, rfftfreq

pd.to_pickle(snakemake, ".fft.py.pkl")
# snakemake = pd.read_pickle(".fft.py.pkl")

with open(f"{snakemake.scriptdir}/to_mne.py", "r") as file:
    exec(file.read())

sfreq = snakemake.params["sfreq"]

sample_data = np.array(pd.read_csv(snakemake.input["vm"], index_col = 0))
episodes = pd.read_csv(snakemake.input["episodes"])

event_fft_freq = []
fft_data = []

for _, episode in episodes.iterrows():
    channel = episode["Channel"]
    start = episode["EventStart"]
    end = int(episode["EventStart"] + sfreq * 0.5)

    if end > episode["EventEnd"]:
        continue

    y = sample_data[channel, start:end]

    fft_res = np.abs(rfft(y))
    event_fft_freq = rfftfreq(len(y), 1 / sfreq)
    fft_data.append(fft_res)

fft_data = np.array(fft_data)
fft_data = pd.DataFrame(fft_data)
fft_data.columns = event_fft_freq

fft_data.to_csv(snakemake.output["fft"], index = False)
