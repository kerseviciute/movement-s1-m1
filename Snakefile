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
        emg = "output/{project}/{animal_id}/{cell_name}/emg/filter.pkl",
        vm = "output/{project}/{animal_id}/{cell_name}/vm/filter.pkl"
    output:
        correlation = "output/{project}/{animal_id}/{cell_name}/lagged_correlation.csv"
    conda: "env/mne.yml"
    script: "python/correlation_lag.py"

#
# Detect episodes of movement
#
rule movement:
    input:
        data = "output/{project}/{animal_id}/{cell_name}/emg/filter.pkl"
    output:
        episodes = "output/{project}/{animal_id}/{cell_name}/movement_episodes.csv"
    params:
        maxTimeApart = config["detect"]["movement"]["maxTimeApart"],
        minEventLength = config["detect"]["movement"]["minLength"],
        percentile = config["detect"]["movement"]["percentile"]
    conda: "env/mne.yml"
    script: "python/movement.py"

#
# Detect episodes of rest
#
rule rest:
    input:
        data = "output/{project}/{animal_id}/{cell_name}/emg/filter.pkl"
    output:
        episodes = "output/{project}/{animal_id}/{cell_name}/rest_episodes.csv"
    params:
        maxTimeApart = config["detect"]["rest"]["maxTimeApart"],
        minEventLength = config["detect"]["rest"]["minLength"],
        percentile = config["detect"]["rest"]["percentile"]
    conda: "env/mne.yml"
    script: "python/rest.py"

rule action_potential:
    input:
        vm = "output/{project}/{animal_id}/{cell_name}/vm/filter.pkl"
    output:
        action_potentials = "output/{project}/{animal_id}/{cell_name}/action_potentials.csv"
    conda: "env/mne.yml"
    script: "python/action_potential.py"
