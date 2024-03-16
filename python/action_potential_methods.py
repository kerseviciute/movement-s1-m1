def find_onset(differential, start_index, min_value = 1):
    value = differential[start_index]
    index = start_index

    while value > min_value and index > 0:
        index -= 1
        value = differential[index]

    return index


def find_offset(differential, start_index, sfreq = 20_000):
    import numpy as np

    value = differential[start_index]
    index = start_index

    # Rising phase - until the AP peak
    while value > 0 and index < len(differential) - 1:
        index += 1
        value = differential[index]

    # Expand the AP duration by 2 ms after peak
    index += round(sfreq * 0.002)
    index = np.min([index, len(differential) - 1])

    return index


def find_ap(signal, time, threshold = 1.5):
    import numpy as np
    import pandas as pd

    differential = np.diff(signal)

    over_threshold = np.where(differential > threshold)[0]

    if len(over_threshold) == 0:
        return []

    idx_diff = np.where(np.diff(over_threshold) != 1)[0]
    detected = np.append(over_threshold[0], over_threshold[idx_diff + 1])

    aps = []
    for ap in detected:
        onset = find_onset(differential, start_index = ap)
        offset = find_offset(differential, start_index = ap)

        aps.append(pd.DataFrame({
            "EventStart": [onset],
            "EventEnd": [offset],
            "Start": [time[onset]],
            "End": [time[onset]],
            "MaxVm": [np.max(signal[onset:offset])]
        }))

    aps = pd.concat(aps)

    return aps
