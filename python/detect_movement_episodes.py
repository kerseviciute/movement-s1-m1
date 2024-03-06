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

def detect_movement_episodes(signal, threshold, sfreq, min_break = 500):
    import numpy as np
    import pandas as pd

    signal_over_threshold = np.abs(signal) > threshold
    change_indices = np.where(np.diff(signal_over_threshold))[0]

    movement_data = pd.DataFrame({
        "EventStart": np.insert(change_indices + 1, 0, 0),
        "EventEnd": np.append(change_indices + 1, len(signal_over_threshold))
    })

    movement_data["Movement"] = signal_over_threshold[movement_data["EventStart"]]
    movement_data = movement_data[movement_data["Movement"]]
    movement_data["EventLength"] = list(movement_data["EventEnd"] - movement_data["EventStart"])

    movement_data = merge_close_events(movement_data, min_break = min_break)

    movement_data["Start"] = movement_data["EventStart"] / sfreq
    movement_data["End"] = movement_data["EventEnd"] / sfreq
    movement_data["Length"] = movement_data["EventLength"] / sfreq

    return movement_data

def detect_rest_episodes(signal, threshold, sfreq, min_break):
    import numpy as np
    import pandas as pd

    signal_under_threshold = np.abs(signal) < threshold
    change_indices = np.where(np.diff(signal_under_threshold))[0]

    movement_data = pd.DataFrame({
        "EventStart": np.insert(change_indices + 1, 0, 0),
        "EventEnd": np.append(change_indices + 1, len(signal_under_threshold))
    })

    movement_data["Movement"] = signal_under_threshold[movement_data["EventStart"]]
    movement_data = movement_data[movement_data["Movement"]]
    movement_data["EventLength"] = list(movement_data["EventEnd"] - movement_data["EventStart"])

    movement_data = merge_close_events(movement_data, min_break = min_break)

    movement_data["Start"] = movement_data["EventStart"] / sfreq
    movement_data["End"] = movement_data["EventEnd"] / sfreq
    movement_data["Length"] = movement_data["EventLength"] / sfreq

    return movement_data
