def merge_close_events(dt, min_break):
    import pandas as pd

    dt = dt.sort_values("EventStart")
    dt = dt.reset_index(drop = True)
    dt["NextStart"] = dt["EventStart"].shift(-1)
    dt["StepBetweenEvents"] = list(dt["NextStart"] - dt["EventStart"])

    i = 0
    while i < dt.shape[0] - 1:
        if dt.iloc[i]["StepBetweenEvents"] >= min_break:
            i += 1
            continue

        current_element = dt.iloc[i]
        next_element = dt.iloc[i + 1]

        # Merge the elements
        dt.at[i, "EventStart"] = current_element["EventStart"]
        dt.at[i, "EventEnd"] = next_element["EventEnd"]

        # Update other parameters
        dt.at[i , "EventLength"] = next_element["EventEnd"] - current_element["EventStart"]
        dt.at[i, "StepBetweenEvents"] = next_element["StepBetweenEvents"]
        dt.at[i, "NextStart"] = next_element["NextStart"]

        dt = dt.drop(i + 1)
        dt = dt.reset_index(drop = True)

    return dt[["EventStart", "EventEnd", "Movement", "EventLength"]]
