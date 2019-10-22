"""
Purpose: Store various simulation parameters in one file.
"""

# standard
import logging

# third party
import numpy as np

# first party
from src.utils.geo.locations import Locations

# seed
SEED = 20190927

# sensors and geography
FIRST_EPIWEEK = 201045
N_TRAIN_WEEKS = 156  # (3 years)
MAX_CACHE_EPIWEEK = 201850

# -- [public] toy simulation set-up --
SENSORS = ["sar3", "epic"]
ATOM_LIST = ["pa", "va", "nc", "sc", "wv"]
REGION_LIST = ATOM_LIST + ["hhs3", "hhs4"] + ["nat"]
EXCLUDE_LOC = list(set(Locations.region_list) - set(REGION_LIST))
# -- [public] end --

# -- [private] assumes access to all data sensors and sources --
# SENSORS = ['gft', 'ght', 'twtr', 'wiki', 'cdc', 'epic', 'sar3', 'arch']
# EXCLUDE_LOC = ["pr", "vi"] + Locations.cen_list
# ATOM_LIST = [l for l in Locations.atom_list if l not in EXCLUDE_LOC]
# REGION_LIST = [l for l in Locations.region_list if l not in EXCLUDE_LOC]
# -- [private] end --

# cross-validation
N_CV_TIMEPOINTS = 10
RIDGE_PARAMS = list(np.exp(np.linspace(np.log(10), np.log(300), 20)))
RIDGE_PARAMS.extend([1e-6, 350., 400.])  # add extreme values for completeness
RIDGE_PARAMS = list(sorted(RIDGE_PARAMS))
LASSO_PARAMS = list(np.exp(np.linspace(0, np.log(2.5), 15)) - 1 + 1e-6)
LASSO_PARAMS.extend([0.001, 0.005, 0.01, 2.])  # add extreme values
LASSO_PARAMS = list(sorted(LASSO_PARAMS))

# random forest parameters
N_ESTIMATORS = 200

# logging
logging.basicConfig(level=logging.INFO)
