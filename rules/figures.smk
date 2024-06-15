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
