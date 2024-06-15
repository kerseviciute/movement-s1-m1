#
# Convert the raw VM .txt files to MNE objects for further processing.
#
rule vm_csv:
    input:
        raw = "raw/{sid}_{cell}_Vm_rest.txt"
    output:
        raw = temp("output/{project}/{sid}/{cell}/vm/raw.csv")
    conda: "../env/mne.yml"
    script: "../python/create_csv.py"

#
# Filter the VM data.
#
rule vm_filter:
    input:
        raw = "output/{project}/{sid}/{cell}/vm/raw.csv"
    output:
        filter = "output/{project}/{sid}/{cell}/vm/filter.csv"
    params:
        drop_below = config["filter"]["vm"]["drop_below"],
        drop_above = config["filter"]["vm"]["drop_above"],
        scale = config["filter"]["vm"]["scale"],
        ch_type = "bio",
        sampling_rate = config["sampling_rate"]
    conda: "../env/mne.yml"
    script: "../python/filter.py"

#
# Detect episodes of action potentials
#
rule action_potential:
    input:
        vm = "output/{project}/{animal_id}/{cell_name}/vm/filter.csv"
    output:
        action_potentials = "output/{project}/{animal_id}/{cell_name}/action_potentials.csv"
    params:
        diffThreshold = config["detect"]["ap"]["diffThreshold"],
        minReachedVoltage = config["detect"]["ap"]["minReachedVoltage"],
        sfreq = config["sampling_rate"]
    conda: "../env/mne.yml"
    script: "../python/action_potential.py"
