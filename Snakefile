import pandas as pd

configfile: "config.yml"

samples = pd.read_csv(config["sample_sheet"])

include: "rules/emg.smk"
include: "rules/vm.smk"
include: "rules/reports.smk"

rule all:
    input:
        expand("{deploy_directory}/{page}.html",
            deploy_directory = config["deploy_directory"],
            page = config["report"]["pages"]
        )

#
# Generate a sample sheet from the Excel data.
#
rule sample_sheet:
    input:
        data = "raw/cells patched_naive.xlsx"
    output:
        sample_sheet = config["sample_sheet"],
        samples = expand("output/{project}/samples.RDS", project = config["project"])
    conda: "env/r.yml"
    script: "R/sample_sheet.R"

#
# Lagged correlation between EMG and membrane potential
#
rule correlation:
    input:
        emg = "output/{project}/{animal_id}/{cell_name}/emg/filter.csv",
        vm = "output/{project}/{animal_id}/{cell_name}/vm/filter.csv"
    output:
        correlation = "output/{project}/{animal_id}/{cell_name}/lagged_correlation.csv"
    params:
        sfreq = config["sampling_rate"]
    conda: "env/mne.yml"
    script: "python/correlation_lag.py"

#
# Basic Vm statistics
#
rule vm_statistics:
    input:
        vm = "output/{project}/{animal_id}/{cell_name}/vm/filter.csv",
        action_potentials = "output/{project}/{animal_id}/{cell_name}/action_potentials.csv",
        movement = "output/{project}/{animal_id}/{cell_name}/movement_episodes.csv",
        rest = "output/{project}/{animal_id}/{cell_name}/rest_episodes.csv"
    output:
        statistics = "output/{project}/{animal_id}/{cell_name}/vm_statistics.csv"
    conda: "env/mne.yml"
    script: "python/vm_statistics.py"
