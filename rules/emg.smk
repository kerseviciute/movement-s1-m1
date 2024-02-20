#
# Convert the raw EMG .txt files to MNE objects for further processing.
#
rule emg_mne:
    input:
        raw = "raw/{sid}_{cell}_EMG_rest.txt"
    output:
        raw = "output/{project}/{sid}/{cell}/emg/raw.pkl"
    params:
        sampling_rate = config["sampling_rate"],
        ch_type = "emg"
    conda: "../env/mne.yml"
    script: "../python/create_mne.py"

#
# Filter the EMG data.
#
rule emg_filter:
    input:
        raw = "output/{project}/{sid}/{cell}/emg/raw.pkl"
    output:
        filter = "output/{project}/{sid}/{cell}/emg/filter.pkl"
    params:
        drop_below = config["filter"]["emg"]["drop_below"],
        drop_above = config["filter"]["emg"]["drop_above"],
        scale = config["filter"]["emg"]["scale"],
        ch_type = "emg"
    conda: "../env/mne.yml"
    script: "../python/filter.py"
