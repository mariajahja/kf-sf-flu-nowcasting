"""
Purpose: Store helper functions for simulation in neurips_main.py.
"""
import itertools
import pickle

from src.config import *
from src.utils.epiweek import range_epiweeks


class ResultFile:
    """Store results from simulation."""

    def __init__(self, filename):
        self.filename = filename

    def save_prediction(self, prediction):
        """Pickle prediction, note this appends results."""
        with open(self.filename, 'ab') as f:
            pickle.dump(prediction, f, pickle.HIGHEST_PROTOCOL)
        f.close()


def mean_impute(a):
    """Impute nan with column mean."""
    col_mean = np.nanmean(a, axis=0)
    inds = np.where(np.isnan(a))
    a[inds] = np.take(col_mean, inds[1])
    return a


def cache(ds, sensors=SENSORS, max_ew=MAX_CACHE_EPIWEEK):
    """Cache sensor and wILI values beforehand.

    Note this fits SF on finalized wILI. In practice, wILI values are often
    revised and this caching procedure is not possible.
    """
    inputs = list(itertools.product(sensors, REGION_LIST))
    cache_weeks = list(range_epiweeks(FIRST_EPIWEEK, max_ew, inclusive=False))
    for col, (sen, loc) in enumerate(inputs):
        ds.get_sensor_values(tuple(cache_weeks), loc, sen)
        ds.get_truth_values(tuple(cache_weeks), loc)
        if col % 50 == 0: logging.info(f"Finished {col}/{len(inputs)}")
