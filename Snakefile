import pandas as pd

configfile: "config.yml"

include: "rules/emg.smk"
include: "rules/vm.smk"

samples = pd.read_csv(config["sampleSheet"])

rule all:
    input:
        expand(expand("output/{project}/{sid}/emg/raw.pkl",
            project = config["project"],
            sid = samples["Location"])
        )

#
# Generate a sample sheet from the Excel data.
#
rule sample_sheet:
    input:
        data = "raw/cells patched_naive.xlsx"
    output:
        sample_sheet = config["sampleSheet"],
        samples = expand("output/{project}/samples.RDS", project = config["project"])
    conda: "env/r.yml"
    script: "R/sample_sheet.R"
