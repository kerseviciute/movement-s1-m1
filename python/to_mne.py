def read_csv_as_mne(filename, sfreq, ch_type = "emg"):
    import mne
    import pandas as pd

    data = pd.read_csv(filename, index_col = 0)
    trials = list(data.index.values)

    ch_types = [ch_type for _ in range(0, len(trials))]
    info = mne.create_info(ch_names = trials, sfreq = sfreq, ch_types = ch_types)

    return mne.io.RawArray(data, info)
