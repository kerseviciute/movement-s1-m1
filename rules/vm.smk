#
# Convert the raw VM .txt files to MNE objects for further processing.
#
rule vm_mne:
    input:
        raw = "raw/{sid}_{cell}_Vm_rest.txt"
    output:
        raw = "output/{project}/{sid}/{cell}/vm/raw.pkl"
    params:
        samplingRate = config["samplingRate"],
        ch_type = "bio"
    conda: "../env/mne.yml"
    script: "../python/create_mne.py"
