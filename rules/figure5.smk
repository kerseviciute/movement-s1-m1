
rule figure5_dataB:
    input:
        samples = config["sample_sheet"],
        correlation = expand("output/{project}/{sid}/lagged_correlation.csv",
            project = config["project"],
            sid = samples["Location"]
        )
    output:
        data = "output/{project}/figure5/fig5B_S1_L23_EMG_VM_cross_correlation.csv"
    conda: "../env/r.yml"
    script: "../R/figure5/dataB.R"
