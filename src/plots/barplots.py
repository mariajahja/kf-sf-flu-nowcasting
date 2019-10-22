"""Create barplot figure of results by season."""

# standard
import pickle
from pathlib import Path

# third party
from matplotlib import pyplot as plt

# first party
from src.config import *
from src.utils.delphi_epidata import Epidata
from src.utils.flu_data_source import FluDataSource
from src.utils.geo.locations import Locations

# path to results directory produced from main simulation
# results are produced by flu season (script loops through entire directory)
RESULT_PATH = Path.cwd() / "results"

if __name__ == "__main__":
    ds = FluDataSource(Epidata, [], REGION_LIST)
    ds.signal_key = 'wili'
    ds.cache_key = 'ilinet'
    files = list(sorted(RESULT_PATH.glob('*.p')))
    n_years = len(files)

    # gather errors (looping by year)
    errors = {}
    for j, file in enumerate(files):
        print(file)
        data = []
        with open(file, 'rb') as fr:
            try:
                while True:
                    data.append(pickle.load(fr))
            except EOFError:
                pass
        fr.close()

        year = file.stem
        method_keys = list(sorted(data[0]['preds'][data[0]['ew']].keys()))
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

        groups = [
            ("States",
             [i for i, loc in enumerate(locs) if loc in Locations.atom_list]),
            ("Regions",
             [i for i, loc in enumerate(locs) if loc in Locations.hhs_list]),
            ("National",
             [i for i, loc in enumerate(locs) if loc in Locations.nat_list])
        ]

        for k, (level, idxs) in enumerate(groups):
            if level not in errors:
                errors[level] = {}

            for i, method in enumerate(method_keys):
                if method not in errors[level]:
                    errors[level][method] = {"mae": np.empty(n_years),
                                             "mad": np.empty(n_years)}

                error = results[method][:, idxs] - wili[:, idxs]
                t, _ = error.shape
                errors[level][method]["mae"][j] = np.mean(np.abs(error))
                errors[level][method]["mad"][j] = np.mean(
                    np.abs(error - np.mean(error))) / np.sqrt(t)

    # produce barplot
    plt.rc('font', size=12)  # controls default text sizes
    plt.rc('ytick', labelsize=11)  # fontsize of the tick labels
    barWidth = 0.08

    level = "States"  # "Regions", "National"
    titles = [f"20{str(f.stem)[:2]}-{str(f.stem[2:4])}" for f in files]
    methods_to_plot = ["sf_l2", "sf_l1", "ridge", "lasso",
                       "sf", "rf_sensor", "rf_source"]
    methods_legend = ["SF + L2", "SF + L1",
                      "Ridge", "Lasso", "SF", "RF (sensors)", "RF (sources)"]
    colors = ["#1da161", "#cf6e2d", "#635cbd", "#cf46aa", "#75b009",
              "#eda13e", "#444b69", "#785722"]
    darker_colors = ["#126b40", "#6b3815", "#302d5c", "#782261", "#416107",
                     "#7d4e10", "#1a1d29", "#362710"]

    plt.figure(figsize=(14, 4))
    x_idx = np.arange(n_years)
    start_idx = []
    start_height = []
    for i, method in enumerate(methods_to_plot):
        mae = errors[level][method]["mae"]
        mad = np.add(errors[level][method]["mae"], errors[level][method]["mad"])
        plt.bar(x_idx, mae, width=barWidth, edgecolor="white",
                label=methods_legend[i], color=colors[i],
                yerr=errors[level][method]["mad"], capsize=2,
                ecolor=darker_colors[i])
        plt.bar(x_idx, mad, edgecolor='white',
                width=barWidth, color=colors[i], alpha=0.2)

        start_idx.append(x_idx[0])
        start_height.append(mad[0])
        x_idx = [x + barWidth for x in x_idx]

    for x, h, l in zip(start_idx, start_height, methods_legend):
        if h != h: h = 0
        # magic numbers to center text
        plt.text(x + (barWidth / 1.9), h + 0.05, l, size="smaller",
                 ha='right', va='bottom', rotation=90)

    plt.xlabel('Season')
    plt.ylabel("Mean Absolute Error")
    plt.xticks([r + (3 * barWidth) for r in range(n_years)], titles)
    plt.title(level, size="medium")
    plt.savefig(RESULT_PATH / 'states.pdf', bbox_inches='tight')
    plt.close()
    plt.clf()
