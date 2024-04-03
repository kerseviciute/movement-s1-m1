#
# Convert the raw EMG .txt files to MNE objects for further processing.
#
rule emg_csv:
    input:
        raw = "raw/{sid}_{cell}_EMG_rest.txt"
    output:
        raw = "output/{project}/{sid}/{cell}/emg/raw.csv"
    conda: "../env/mne.yml"
    script: "../python/create_csv.py"

#
# Filter the EMG data.
#
rule emg_filter:
    input:
        raw = "output/{project}/{sid}/{cell}/emg/raw.csv"
    output:
        filter = "output/{project}/{sid}/{cell}/emg/filter.csv"
    params:
        drop_below = config["filter"]["emg"]["drop_below"],
        drop_above = config["filter"]["emg"]["drop_above"],
        scale = config["filter"]["emg"]["scale"],
        ch_type = "emg",
        sampling_rate = config["sampling_rate"]
    conda: "../env/mne.yml"
    script: "../python/filter.py"

#
# Detect episodes of movement
#
rule movement:
    input:
        data = "output/{project}/{animal_id}/{cell_name}/emg/filter.csv"
    output:
        episodes = "output/{project}/{animal_id}/{cell_name}/movement_episodes.csv"
    params:
        maxTimeApart = config["detect"]["movement"]["maxTimeApart"],
        minEventLength = config["detect"]["movement"]["minLength"],
        percentile = config["detect"]["movement"]["percentile"],
        sfreq = config["sampling_rate"]
    conda: "../env/mne.yml"
    script: "../python/movement.py"

#
# Detect episodes of rest
#
rule rest:
    input:
        data = "output/{project}/{animal_id}/{cell_name}/emg/filter.csv"
    output:
        episodes = "output/{project}/{animal_id}/{cell_name}/rest_episodes.csv"
    params:
        maxTimeApart = config["detect"]["rest"]["maxTimeApart"],
        minEventLength = config["detect"]["rest"]["minLength"],
        percentile = config["detect"]["rest"]["percentile"],
        sfreq = config["sampling_rate"]
    conda: "../env/mne.yml"
    script: "../python/rest.py"
