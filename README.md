# Movement in S1 and M1 cortices

Click [here](https://kerseviciute.github.io/movement-s1-m1/) to see full analysis.

## How to run this analysis

The data should be available under the `raw` directory in this repository root (create a symbolic link).

1. Install conda and mamba. Create the _snakemake_ (`env/snakemake.yml`) conda environment
   (required to run the analysis) and the _spontaneous-movement-mne_ (`env/mne.yml`) conda environment
   (required to generate the reports).

```shell
# In macOS M2
export CONDA_SUBDIR=osx-arm64
```

```shell
mamba env create -f env/snakemake.yml
mamba env create -f env/mne.yml
```

2. Activate the snakemake environment.

```shell
conda activate snakemake
```

3. Run ``snakemake``:

```shell
snakemake --conda-frontend mamba --use-conda --rerun-triggers mtime --cores 1 -p all
```

## TODO

### Feb 21

- [x] correlation peak (when? / per cell)
- [x] statistical test on 0 / max timepoint (mean per cell) (compare different groups, e.g. L2/3 vs L5, S1 vs M1)
- [x] emg detection with 10th percentile + compare with TKEO
- [x] high pass filter (2 Hz) - average abs EMG on detected events

### Mar 06

Use movements detected with 10th upper percentile method. Do not filter out any
detected episodes for now.

- [x] show that correlation is different from zero
- [x] color by cortex in correlation plots
- [x] detect EMG off
  - [ ] What is happening in W3? How to filter out the heartbeat?
- [ ] ap detection (90 percentile?)
- [ ] compare movement on/off (mean vm, sd, ap)
  - [ ] frequency analysis
- [ ] kursinio planas

### Mar ?

- [ ] compare conductance during on/off
- [ ] EMG event filtering
