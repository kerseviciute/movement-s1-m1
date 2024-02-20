#
# Convert the raw EMG .txt files to MNE objects for further processing.
#
rule emg_mne:
    input:
        raw = "raw/{sid}_{cell}_EMG_rest.txt"
    output:
        raw = "output/{project}/{sid}/{cell}/emg/raw.pkl"
    params:
        samplingRate = config["samplingRate"],
        ch_type = "emg"
    conda: "../env/mne.yml"
    script: "../python/create_mne.py"
