"""
Purpose: Run nowcasting simulation for methods listed in the paper:

Kalman Filter, Sensor Fusion, and Constrained Regression: Equivalences
and Insights
"""

# standard
import datetime
from multiprocessing import Pool, cpu_count

# third party
import click
from sklearn.ensemble import RandomForestRegressor

# first party
from src.models import sf, reg
from src.utils.delphi_epidata import Epidata
from src.utils.epiweek import add_epiweeks
from src.utils.flu_data_source import FluDataSource
from src.utils.sim_helper import *
from src.utils.us_fusion import UsFusion


def get_training_data(ew, ds, n_train_weeks=N_TRAIN_WEEKS, regions=REGION_LIST):
    inputs = list(itertools.product(SENSORS, regions))
    n_inputs = len(inputs)

    test_weeks = list([ew])
    train_weeks = list(range_epiweeks(FIRST_EPIWEEK, ew, inclusive=False))
    if len(train_weeks) >= n_train_weeks:
        train_weeks = train_weeks[len(train_weeks) - n_train_weeks:]
    n_train_weeks = len(train_weeks)
    logging.info(f"Number of training weeks is {n_train_weeks}")

    readings = np.full((1, n_inputs), np.nan)
    sensor_vals = np.full((n_train_weeks, n_inputs), np.nan)
    for col, (sen, loc) in enumerate(inputs):
        sensor_vals[:, col] = ds.get_sensor_values(tuple(train_weeks), loc, sen)

        for row, week in enumerate(test_weeks):
            value = ds.get_sensor_value(week, loc, sen)
            if value is not None:
                readings[row, col] = value

    hist_wili = np.full((n_train_weeks, len(ATOM_LIST)), np.nan)
    for col, loc in enumerate(ATOM_LIST):
        hist_wili[:, col] = ds.get_truth_values(tuple(train_weeks), loc)

    # remove empty columns
    get_finite_columns = lambda d: np.any(np.isfinite(d), axis=0)
    finite_sensors = get_finite_columns(sensor_vals)
    finite_readings = get_finite_columns(readings)
    keep_columns = np.logical_and(finite_sensors, finite_readings)
    inputs = list(itertools.compress(inputs, keep_columns))
    sensor_vals = sensor_vals[:, keep_columns]
    readings = readings[:, keep_columns]

    # construct H, W
    selected_inputs = tuple([loc for (sens, loc) in inputs])
    H, W, output_locs = UsFusion.determine_statespace(
        selected_inputs, season=None, exclude_locations=tuple(EXCLUDE_LOC))

    # if any current readings are nan, then fill with historical mean
    if np.sum(np.isnan(readings)) != 0:
        means = np.nanmean(sensor_vals, axis=0)
        inds = np.where(np.isnan(readings))
        readings[inds] = np.take(means, inds[1])

    logging.debug(hist_wili.shape)
    logging.info(f"Shape of Z is {sensor_vals.shape}")
    return {"wili": mean_impute(hist_wili),
            "sensors": mean_impute(sensor_vals),
            "new_sensors": readings,
            "H": H, "W": W, "output_locs": output_locs}


def cv(method_key, cv_ew, method, X, Z, H, new_z, truth, params):
    """Cross-validation function for parallelization."""
    logging.debug(f'[CV] Running {method_key} for {cv_ew}')
    Beta = method(X, Z, H, params)
    x_hats = new_z @ Beta
    errors = np.mean(np.abs(np.subtract(x_hats, truth)), axis=1)
    return {"cv_ew": cv_ew, "method_key": method_key, "errors": errors}


def predict(method_key, ew, method, X, Z, H, W, new_z, param):
    """Prediction function for parallelization."""
    logging.debug(f'[FINAL] Running {method_key} for {ew}'
                  f'. The best param is {param:.3f}.')
    Beta = method(X, Z, H, [param])
    x_hat = ((new_z @ Beta) @ W.T).flatten()
    return {"x_hat": x_hat, "method_key": method_key}


def run(ew_to_pred, cv_dict, methods, ds):
    data = get_training_data(ew_to_pred, ds)
    assert np.sum(np.isnan(data["sensors"])) == 0

    # calculate one-week-ahead prediction error for parameter grid
    # this is like cross-validation, for this time series context
    cv_weeks = list(range_epiweeks(add_epiweeks(ew_to_pred, -N_CV_TIMEPOINTS),
                                   ew_to_pred, inclusive=False))
    for i, cv_ew in enumerate(cv_weeks):
        if cv_ew in cv_dict:  # if stored, skip calculation
            logging.debug(f"Found stored result for {cv_ew}, skipping")
            continue

        end_week_idx = (data["wili"].shape[0] - (i + 1))
        X = data["wili"][:end_week_idx, :]
        Z = data["sensors"][:end_week_idx, :]
        H = data["H"]
        new_z = data["sensors"][end_week_idx, :]
        truth = data["wili"][end_week_idx, :]

        pool_results = []
        for method_key, method, params in methods:
            pool_results.append(pool.apply_async(cv, args=(
                method_key, cv_ew, method, X, Z, H, new_z, truth, params)))

        pool_results = [proc.get() for proc in pool_results]
        for res in pool_results:
            cv_ew = res["cv_ew"]
            if cv_ew not in cv_dict: cv_dict[cv_ew] = {}
            cv_dict[cv_ew][res["method_key"]] = res["errors"]

    # start prediction
    predictions = {ew_to_pred: {}}

    # run regularized methods with chosen cv parameters
    pool_results = []
    for method_key, method, params in methods:
        errors = np.array([cv_dict[ew][method_key] for ew in cv_weeks])
        best_param = params[np.argmin(np.mean(errors, axis=0))]

        pool_results.append(pool.apply_async(predict, args=(
            method_key, ew_to_pred, method, data["wili"], data["sensors"],
            data["H"], data["W"], data["new_sensors"], best_param)))

    pool_results = [proc.get() for proc in pool_results]
    for res in pool_results:
        predictions[ew_to_pred][res["method_key"]] = res["x_hat"]

    # run sf with no regularization
    logging.info(f"[FINAL] Running sf for {ew_to_pred}.")
    sf_Beta = sf.sf_l2(data["wili"], data["sensors"], data["H"], [0])
    sf_x_hat = ((data["new_sensors"] @ sf_Beta) @ data["W"].T).flatten()
    predictions[ew_to_pred]["sf"] = sf_x_hat

    # run regression with no regularization
    logging.info(f"[FINAL] Running regression for {ew_to_pred}.")
    reg_Beta = reg.ridge(data["wili"], data["sensors"], data["H"], [0])
    reg_x_hat = ((data["new_sensors"] @ reg_Beta) @ data["W"].T).flatten()
    predictions[ew_to_pred]["reg"] = reg_x_hat

    # run random forest (sensors)
    logging.info(f"[FINAL] Running RF for {ew_to_pred}.")
    max_features = int(data["sensors"].shape[1] / 3)
    all_results = np.empty(len(ATOM_LIST))
    for i in range(len(ATOM_LIST)):
        rfc = RandomForestRegressor(n_estimators=N_ESTIMATORS,
                                    max_features=max_features)
        rfc.fit(data["sensors"], data["wili"][:, i])
        all_results[i] = rfc.predict(data["new_sensors"].reshape(1, -1))
    rf_x_hat = (data["W"] @ all_results.reshape(-1, 1)).flatten()
    predictions[ew_to_pred]["rf_sensor"] = rf_x_hat

    return {"ew": ew_to_pred, "preds": predictions, "locs": data["output_locs"]}


@click.command()
@click.argument('start', type=int)
@click.argument('end', type=int)
@click.argument('out', type=str)
def init(start, end, out):
    inputs = list(itertools.product(SENSORS, REGION_LIST))
    ds = FluDataSource(Epidata, SENSORS, inputs)  # FluDataSource on Delphi side
    ds.signal_key = 'wili'
    ds.cache_key = 'ilinet'
    cache(ds)  # cache sensors for efficiency

    cv_methods = [("sf_l2", sf.sf_l2, RIDGE_PARAMS),
                  ("sf_l1", sf.sf_l1, LASSO_PARAMS),
                  ("ridge", reg.ridge, RIDGE_PARAMS),
                  ("lasso", reg.lasso, LASSO_PARAMS)]
    cv_dict = {}

    filename = datetime.datetime.now().strftime("%Y%m%d.p")
    out_file = ResultFile(out + "-" + filename)

    for ew in [ew for ew in range_epiweeks(start, end)]:
        pred = run(ew, cv_dict, cv_methods, ds)
        logging.debug(pred)
        out_file.save_prediction(pred)
        logging.info(f"Finished {ew}.")


if __name__ == "__main__":
    np.random.seed(SEED)

    # set-up multiprocessing
    n_cpu = min(4, cpu_count())  # going through 4 methods for cv
    pool = Pool(n_cpu)
    logging.info(f'Starting job, using {n_cpu} CPUs.')

    init()
