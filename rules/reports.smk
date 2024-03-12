rule report_summary:
    output:
        report = "{deploy_directory}/index.html"
    params:
        script = "reports/index.Rmd",
    conda: "../env/r.yml"
    script: "../R/render.R"

rule report_correlation:
    input:
        samples = config["sample_sheet"],
        correlation = expand("output/{project}/{sid}/lagged_correlation.csv",
            project = config["project"],
            sid = samples["Location"]
        )
    output:
        report = "{deploy_directory}/correlation.html"
    params:
        script = "reports/correlation.Rmd",
        prefix = expand("output/{project}", project = config["project"])
    conda: "../env/r.yml"
    script: "../R/render.R"

rule report_preprocess:
    input:
        samples = config["sample_sheet"],
        emg_raw = expand("output/{project}/{sid}/emg/raw.pkl",
            project = config["project"],
            sid = samples["Location"]),
        emg_filter = expand("output/{project}/{sid}/emg/filter.pkl",
            project = config["project"],
            sid = samples["Location"]),
        vm_filter = expand("output/{project}/{sid}/vm/filter.pkl",
            project = config["project"],
            sid = samples["Location"])
    output:
        report = "{deploy_directory}/preprocess_{region}.html"
    params:
        script = "reports/preprocess.Rmd",
        prefix = expand("output/{project}", project = config["project"])
    conda: "../env/r.yml"
    script: "../R/render.R"

rule report_emg:
    input:
        samples = config["sample_sheet"],
        emg_raw = expand("output/{project}/{sid}/emg/raw.pkl",
            project = config["project"],
            sid = samples["Location"]),
        emg_filter = expand("output/{project}/{sid}/emg/filter.pkl",
            project = config["project"],
            sid = samples["Location"]),
        movement = expand("output/{project}/{sid}/movement_episodes.csv",
            project = config["project"],
            sid = samples["Location"]),
        rest = expand("output/{project}/{sid}/rest_episodes.csv",
            project = config["project"],
            sid = samples["Location"])
    output:
        report = "{deploy_directory}/emg_{region}.html"
    params:
        script = "reports/emg.Rmd",
        prefix = expand("output/{project}", project = config["project"])
    conda: "../env/r.yml"
    script: "../R/render.R"
