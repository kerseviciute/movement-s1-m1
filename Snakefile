import pandas as pd

configfile: "config.yml"

include: "rules/emg.smk"
include: "rules/vm.smk"

samples = pd.read_csv(config["sample_sheet"])

rule all:
    input:
        expand("output/{project}/{sid}/lagged_correlation.csv",
            project = config["project"],
            sid = samples["Location"]
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
        emg = "output/{{project}}/{animal_id}/{cell_name}/emg/filter.pkl",
        vm = "output/{{project}}/{animal_id}/{cell_name}/vm/filter.pkl"
    output:
        correlation = "output/{project}/{animal_id}/{cell_name}/lagged_correlation.csv"
    conda: "env/mne.yml"
    script: "python/correlation.py"
