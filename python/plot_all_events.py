def plot_all_events(data, channels, movement = None, no_movement = None, limit = None, show_full = True, alpha = 1):
    from matplotlib import pyplot as plt
    import numpy as np

    if limit is None:
        limit = [-4, 4]

    limit_start = limit[0]
    limit_end = limit[1]

    figure, axes = plt.subplots(nrows = len(data.ch_names), figsize = (20, 0.5 * len(data.ch_names)))
    plt.subplots_adjust(
        left = 0.1,
        bottom = 0.1,
        right = 0.9,
        top = 0.9,
        wspace = 0,
        hspace = 0
    )

    x = data.times

    for i, channel in enumerate(channels):
        y = data.get_data(picks = [i])[0]

        if show_full:
            axes[i].plot(x, y, linewidth = 0.5, color = "black")

        if limit is not None:
            axes[i].set_ylim(limit_start, limit_end)

        axes[i].get_xaxis().set_ticks([])
        axes[i].get_yaxis().set_ticks([])
        axes[i].set_ylabel(channel, rotation = 0, labelpad = 60, loc = "center")
        axes[i].set_xlim(0, np.max(data.times))

        if movement is not None:
            events = movement[movement["Channel"] == i]
            for index, row in events.iterrows():
                start = int(row["EventStart"])
                end = int(row["EventEnd"])

                axes[i].plot(x[start:end], y[start:end], linewidth = 1, color = "red", alpha = alpha)

        if no_movement is not None:
            events = no_movement[no_movement["Channel"] == i]
            for index, row in events.iterrows():
                start = int(row["EventStart"])
                end = int(row["EventEnd"])

                axes[i].plot(x[start:end], y[start:end], linewidth = 1, color = "blue", alpha = alpha)

        if i != 0:
            axes[i].spines["top"].set_visible(False)

        if i != len(channels) - 1:
            axes[i].spines["bottom"].set_visible(False)

    plt.show()
