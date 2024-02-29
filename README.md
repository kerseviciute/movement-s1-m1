# Movement in S1 and M1 cortices

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
snakemake --conda-frontend mamba --use-conda --cores 1 -p all
```

## TODO

### Feb 21

- [x] correlation peak (when? / per cell)
- [x] statistical test on 0 / max timepoint (mean per cell) (compare different groups, e.g. L2/3 vs L5, S1 vs M1)
- [ ] emg detection with 10th percentile + compare with TKEO
- [ ] high pass filter (2 Hz) - average abs EMG on detected events
