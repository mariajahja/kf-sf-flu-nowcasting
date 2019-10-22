"""Create prediction figure with underlying ILI level."""

# standard
import datetime
import pickle
from pathlib import Path

# third party
from matplotlib import pyplot as plt
from matplotlib.dates import DateFormatter
from matplotlib.ticker import FormatStrFormatter

# first party
from src.config import *
from src.utils.delphi_epidata import Epidata
from src.utils.flu_data_source import FluDataSource

# path to results file produced from main simulation
# results are produced by flu season
RESULT_FILE_PATH = Path.cwd() / "results" / "1718.p"

if __name__ == "__main__":
    ds = FluDataSource(Epidata, SENSORS, REGION_LIST)
    ds.signal_key = 'wili'
    ds.cache_key = 'ilinet'

    # get nowcasts
    data = []
    with open(RESULT_FILE_PATH, 'rb') as fr:
        try:
            while True:
                data.append(pickle.load(fr))
        except EOFError:
            pass

    method_keys = data[0]['preds'][data[0]['ew']].keys()
    locs = data[0]['locs']
    n_weeks = len(data)
    n_locs = len(locs)

    results = {}
    weeks = []
    for i, d in enumerate(data):
        ew = d['ew']
        weeks.append(ew)
        assert set(d['preds'][ew].keys()) == set(method_keys)
        assert d['locs'] == locs

        for method in method_keys:
            if method not in results:
                results[method] = np.empty((n_weeks, n_locs))
            results[method][i, :] = d['preds'][ew][method]

    wili = np.empty((n_weeks, n_locs))
    for i, loc in enumerate(locs):
        wili[:, i] = ds.get_truth_values(tuple(weeks), loc)

    # reformat epiweek to date for plotting
    pw = []
    for week in weeks:
        d = str(int(week) + 1)
        pw.append(datetime.datetime.strptime(d + '-5', '%Y%W-%w'))

    ## start plot
    plt.rc('font', size=16)  # controls default text sizes
    plt.rc('ytick', labelsize=14)  # fontsize of the tick labels

    green = "seagreen"
    orange = "darkorange"
    purple = "darkviolet"
    pink = "magenta"
    blue = "royalblue"

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    loc_idx = np.where(np.array(locs) == 'sc')[0][0]
    ax.plot(pw, wili[:, loc_idx], 'black', label="wILI", linewidth=2.5,
            alpha=.9)
    ax.plot(pw, results["sf_l2"][:, loc_idx], pink, marker='*',
            label="SF + L2", linewidth=1.2, alpha=0.7)
    ax.plot(pw, results["ridge"][:, loc_idx], blue, marker='+',
            label="Ridge", linewidth=2, alpha=0.9, linestyle='--')
    ax.plot(pw, results["sf"][:, loc_idx], orange, linestyle='--',
            label="SF", linewidth=1.2, alpha=0.9)
    ax.plot(pw, results["rf_sensor"][:, loc_idx], c=green,
            label="RF (sensors)", marker="d", markersize=4)
    ax.xaxis.set_major_formatter(DateFormatter('%b'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.legend(loc='upper right')

    plt.tight_layout()

    plt.savefig(RESULT_FILE_PATH.parent / "sc_1718.pdf")
    plt.show()
    plt.clf()
