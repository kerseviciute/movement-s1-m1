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
        vm_filter = expand("output/{project}/{sid}/vm/filter.pkl",
            project = config["project"],
            sid = samples["Location"]),
        action_potentials = expand("output/{project}/{sid}/action_potentials.csv",
            project = config["project"],
            sid = samples["Location"])
    output:
        report = "{deploy_directory}/vm_{region}.html"
    params:
        script = "reports/vm.Rmd",
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

rule report_movement_vs_rest:
    input:
        samples = config["sample_sheet"],
        movement = expand("output/{project}/{sid}/movement_episodes.csv",
            project = config["project"],
            sid = samples["Location"]),
        rest = expand("output/{project}/{sid}/rest_episodes.csv",
            project = config["project"],
            sid = samples["Location"]),
        action_potentials = expand("output/{project}/{sid}/action_potentials.csv",
            project = config["project"],
            sid = samples["Location"])
    output:
        report = "{deploy_directory}/movement_vs_rest.html"
    params:
        script = "reports/movement_vs_rest.Rmd",
        prefix = expand("output/{project}", project = config["project"])
    conda: "../env/r.yml"
    script: "../R/render.R"