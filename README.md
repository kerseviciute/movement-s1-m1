# Movement in S1 and M1 cortices

**Complete up-to-date analysis is available [here](https://kerseviciute.github.io/movement-s1-m1/).**

The Jupyter notebooks (under the `notebooks` directory) are used primarily for algorithm and method
development and usually include runs for a single sample. After, the finalized algorithms
are moved to the pipeline (to the `R` or `python` directories) and are used for all samples.
While I try to rerun the notebooks for all changes, sometimes **the notebooks may not be up-to-date**.
Please see the `docs` directory or the [deployed website](https://kerseviciute.github.io/movement-s1-m1/)
for the most recent results.

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
  - [x] What is happening in W3? How to filter out the heartbeat?
  - [x] Add counts of EMG on / off in the report (number of events, total time)
- [x] AP detection (90 percentile?)
  - [x] filter out false positives!
- [x] vm report html
- [x] compare movement on/off
  - [x] number of episodes
  - [x] episode length
  - [x] mean vm
  - [x] vm sd
  - [x] AP
- [x] correlation report html

**To discuss**
- Filtering had to be changed because heart beat was disturbing to the rest period detection. In some cases, it
  still is, and I am not sure how to get rid of it.
  - Few samples have very low number of rest episodes (e.g. W2 C3 S1 L23).
- Decided to use differential analysis for detection of action potentials (90th percentile and any average-like
  methods will likely introduce a lot of false positives).
  - AP detection: the hell happened in W4 C15 (S1 L5)?
- Why does the Vm continuously grow over time in some samples (e.g. W4 C11 S1 L5)? Is this some unwanted effect
  which should be corrected for?
- Combining S1 layers into a single group seems anti-productive since there are significant differences in the
  correlation patterns between the layers L2/3 and L5 of S1 (not seen in M1, though).

### Mar 20

- [ ] fix correlation report (max corr, non absolute for test)
- [ ] use low pass filtered data for detection
- [ ] dying cells
- [ ] frequency analysis (movement on / movement off)
- [ ] mixed effects models for time?
- [ ] kursinio planas
- [ ] add method descriptions to reports
  - [ ] index (how to use the website)
  - [ ] correlation analysis
  - [ ] EMG (filtering + movement detection)
  - [ ] Vm (filtering + AP detection)

### Mar ?

- [ ] compare conductance during on/off
- [ ] EMG event filtering
