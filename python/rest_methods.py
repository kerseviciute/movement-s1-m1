def detect_rest_episodes(signal, threshold, sfreq, min_break):
    import numpy as np
    import pandas as pd

    signal_over_threshold = np.abs(signal) > threshold
    change_indices = np.where(np.diff(signal_over_threshold))[0]

    movement_data = pd.DataFrame({
        "EventStart": np.insert(change_indices + 1, 0, 0),
        "EventEnd": np.append(change_indices + 1, len(signal_over_threshold))
    })

    movement_data["Movement"] = signal_over_threshold[movement_data["EventStart"]]
    movement_data = movement_data[movement_data["Movement"] == False]
    movement_data["EventLength"] = list(movement_data["EventEnd"] - movement_data["EventStart"])

    movement_data = merge_close_events(movement_data, min_break = min_break)

    movement_data["Start"] = movement_data["EventStart"] / sfreq
    movement_data["End"] = movement_data["EventEnd"] / sfreq
    movement_data["Length"] = movement_data["EventLength"] / sfreq

    return movement_data
