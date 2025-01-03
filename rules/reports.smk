rule report_summary:
    input:
        samples = config["sample_sheet"]
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
        script = "reports/correlation.Rmd"
    conda: "../env/r.yml"
    script: "../R/render.R"

rule report_vm:
    input:
        samples = config["sample_sheet"],
        vm_filter = expand("output/{project}/{sid}/vm/filter.csv",
            project = config["project"],
            sid = samples["Location"]),
        action_potentials = expand("output/{project}/{sid}/action_potentials.csv",
            project = config["project"],
            sid = samples["Location"])
    output:
        report = "{deploy_directory}/vm_{region}.html"
    params:
        script = "reports/vm.Rmd"
    conda: "../env/r.yml"
    script: "../R/render.R"

rule report_emg:
    input:
        samples = config["sample_sheet"],
        emg_filter = expand("output/{project}/{sid}/emg/filter.csv",
            project = config["project"],
            sid = samples["Location"]),
        movement = expand("output/{project}/{sid}/movement_{{type}}.csv",
            project = config["project"],
            sid = samples["Location"]),
        rest = expand("output/{project}/{sid}/rest_{{type}}.csv",
            project = config["project"],
            sid = samples["Location"])
    output:
        report = "{deploy_directory}/{type}_emg_{region}.html"
    params:
        script = "reports/emg.Rmd"
    conda: "../env/r.yml"
    script: "../R/render.R"

rule report_movement_vs_rest:
    input:
        samples = config["sample_sheet"],
        movement = expand("output/{project}/{sid}/movement_{{type}}.csv",
            project = config["project"],
            sid = samples["Location"]),
        rest = expand("output/{project}/{sid}/rest_{{type}}.csv",
            project = config["project"],
            sid = samples["Location"]),
        action_potentials = expand("output/{project}/{sid}/action_potentials.csv",
            project = config["project"],
            sid = samples["Location"]),
        statistics = expand("output/{project}/{sid}/vm_statistics.csv",
            project = config["project"],
            sid = samples["Location"])
    output:
        report = "{deploy_directory}/{type}_movement_vs_rest.html"
    params:
        script = "reports/movement_vs_rest.Rmd"
    conda: "../env/r.yml"
    script: "../R/render.R"

rule report_final_emg_episodes:
    input:
        samples = config["sample_sheet"],

        # Data before filtering
        movement = expand("output/{project}/{sid}/movement_episodes.csv",
            project = config["project"],
            sid = samples["Location"]),
        rest = expand("output/{project}/{sid}/rest_episodes.csv",
            project = config["project"],
            sid = samples["Location"]),

        # Movement onset
        movement_filter_onset = expand("output/{project}/{sid}/movement_filtered_onset.csv",
            project = config["project"],
            sid = samples["Location"]),
        rest_filter_onset = expand("output/{project}/{sid}/rest_filtered_onset.csv",
            project = config["project"],
            sid = samples["Location"]),
        onset_emg = expand("output/{project}/{sid}/emg/movement_onset.csv",
            project = config["project"],
            sid = samples["Location"]),
        onset_vm = expand("output/{project}/{sid}/vm/movement_onset.csv",
            project = config["project"],
            sid = samples["Location"]),

        # Movement offset
        movement_filter_offset = expand("output/{project}/{sid}/movement_filtered_offset.csv",
            project = config["project"],
            sid = samples["Location"]),
        rest_filter_offset = expand("output/{project}/{sid}/rest_filtered_offset.csv",
            project = config["project"],
            sid = samples["Location"]),
        offset_emg = expand("output/{project}/{sid}/emg/movement_offset.csv",
            project = config["project"],
            sid = samples["Location"]),
        offset_vm = expand("output/{project}/{sid}/vm/movement_offset.csv",
            project = config["project"],
            sid = samples["Location"]),

        # Action potentials
        action_potentials = expand("output/{project}/{sid}/action_potentials.csv",
            project = config["project"],
            sid = samples["Location"])
    output:
        report = "{deploy_directory}/analysis_vm_dynamics.html"
    params:
        script = "reports/vm_dynamics.Rmd"
    conda: "../env/r.yml"
    script: "../R/render.R"

rule report_fft:
    input:
        samples = config["sample_sheet"],
        movement = expand("output/{project}/{sid}/movement_fft.csv",
            project = config["project"],
            sid = samples["Location"]),
        rest = expand("output/{project}/{sid}/rest_fft.csv",
            project = config["project"],
            sid = samples["Location"])
    output:
        report = "{deploy_directory}/fft.html"
    params:
        script = "reports/fft.Rmd"
    conda: "../env/r.yml"
    script: "../R/render.R"
