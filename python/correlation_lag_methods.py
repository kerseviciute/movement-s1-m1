def filter_emg(data):
    from scipy import signal
    import numpy as np

    butter_filter = signal.butter(
        N = 20,
        Wn = [ 100 ],
        btype = "lowpass",
        output = "sos",
        fs = 20_000
    )

    return signal.sosfiltfilt(butter_filter, np.abs(data))

def correlation_lag(emg, vm, lag_range = None, lag_step = 50):
    from scipy.stats import pearsonr
    from scipy.ndimage import shift
    import pandas as pd

    if lag_range is None:
        lag_range = [-5000, 5000]

    lags = []
    correlations = []

    for i, channel in enumerate(emg.ch_names):
        emg_data = emg.copy().pick(picks = [i]).get_data()[ 0, : ]
        vm_data = vm.copy().pick(picks = [i]).get_data()[ 0, : ]

        filtered = filter_emg(emg_data)

        x = []
        y = []

        for lag in range(lag_range[0], lag_range[1], lag_step):
            x.append(lag / emg.info["sfreq"])
            y.append(pearsonr(vm_data, shift(filtered, lag)).statistic)

        lags = x
        correlations.append(y)

    return [lags, pd.DataFrame(correlations)]
