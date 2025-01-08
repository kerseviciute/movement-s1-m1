
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

rule figure5_dataC1:
    input:
        emg = expand(
            "output/{{project}}/{animalID}/{cellID}/emg/filter.csv",
            animalID = "W1", cellID = "C2"
        ),
        movement = expand(
            "output/{{project}}/{animalID}/{cellID}/movement_episodes.csv",
            animalID = "W1", cellID = "C2"
        ),
        rest = expand(
            "output/{{project}}/{animalID}/{cellID}/rest_episodes.csv",
            animalID = "W1", cellID = "C2"
        )
    params:
        animalID = "W1",
        cellID = "C2"
    output:
        data = "output/{project}/figure5/fig5C_S1_L23_W1C2_EMG_example.csv"
    conda: "../env/r.yml"
    script: "../R/figure5/dataC1.R"

rule figure5_dataC2:
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
        data = "output/{project}/figure5/fig5C_S1_L23_Vm_Mean_SD.csv",
        model = "output/{project}/figure5/fig5C_S1_L23_Vm_Mean_SD_model_fit.csv"
    conda: "../env/r.yml"
    script: "../R/figure5/dataC2.R"

rule figure5_dataD:
    input:
        samples = config["sample_sheet"],
        movement_filter = expand("output/{project}/{sid}/movement_filtered_onset.csv",
            project = config["project"],
            sid = samples[ samples.Region == "S1_L23" ]["Location"]),
        rest_filter = expand("output/{project}/{sid}/rest_filtered_onset.csv",
            project = config["project"],
            sid = samples[ samples.Region == "S1_L23" ]["Location"]),
        emg = expand("output/{project}/{sid}/emg/movement_onset.csv",
            project = config["project"],
            sid = samples[ samples.Region == "S1_L23" ]["Location"]),
        vm = expand("output/{project}/{sid}/vm/movement_onset.csv",
            project = config["project"],
            sid = samples[ samples.Region == "S1_L23" ]["Location"])
    output:
        average = "output/{project}/figure5/fig5D_S1_L23_vm_average.csv",
        average_period = "output/{project}/figure5/fig5D_S1_L23_vm_period_average.csv",
        model = "output/{project}/figure5/fig5D_S1_L23_vm_period_average_model_fit.csv",
        vm_png = "output/{project}/figure5/fig5D_S1_L23_emg_onset_vm.png",
        cellOnsetVm = "output/{project}/figure5/fig5D_S1_L23_movement_onset_vm_per_cell.csv"
    conda: "../env/r.yml"
    script: "../R/figure5/dataD.R"
