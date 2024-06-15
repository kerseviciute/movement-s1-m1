rule figure1:
    input:
        samples = config["sample_sheet"],
        correlation = expand("output/{project}/{sid}/lagged_correlation.csv",
            project = config["project"],
            sid = samples["Location"]
        )
    output:
        figure = "docs/www/figure1.png"
    conda: "../env/r.yml"
    script: "../R/figures/figure1.R"

rule figure2:
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
            sid = samples["Location"]),
        statistics = expand("output/{project}/{sid}/vm_statistics.csv",
            project = config["project"],
            sid = samples["Location"])
    output:
        figure = "docs/www/figure2.png"
    conda: "../env/r.yml"
    script: "../R/figures/figure2.R"

rule figure3:
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
        movement_filter = expand("output/{project}/{sid}/movement_filtered_onset.csv",
            project = config["project"],
            sid = samples["Location"]),
        rest_filter = expand("output/{project}/{sid}/rest_filtered_onset.csv",
            project = config["project"],
            sid = samples["Location"]),
        emg = expand("output/{project}/{sid}/emg/movement_onset.csv",
            project = config["project"],
            sid = samples["Location"]),
        vm = expand("output/{project}/{sid}/vm/movement_onset.csv",
            project = config["project"],
            sid = samples["Location"]),

        # Action potentials
        action_potentials = expand("output/{project}/{sid}/action_potentials.csv",
            project = config["project"],
            sid = samples["Location"])
    output:
        figure = "{deploy_directory}/www/figure3.png"
    conda: "../env/r.yml"
    script: "../R/figures/figure3.R"
