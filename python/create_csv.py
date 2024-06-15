import pandas as pd

pd.to_pickle(snakemake, ".create_csv.py.pkl")
# snakemake = pd.read_pickle(".create_csv.py.pkl")

# Read raw data
data = pd.read_csv(snakemake.input["raw"], sep = "\t", index_col = False)

# Create CSV object (trials as rows, time as columns)
trials = data.columns
trials = [trial.replace("'", "") for trial in trials]

data = data.transpose()
data.index = trials

data.to_csv(snakemake.output["raw"])
