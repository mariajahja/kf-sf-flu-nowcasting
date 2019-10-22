"""
This file was created by the CMU DELPHI Epi-forecasting group.
https://github.com/cmu-delphi/

===============
=== Purpose ===
===============

A wrapper for the Epidata API as used for nowcasting. Caching is used
extensively to reduce the number of requests made to the API.
"""

# standard
import abc
import functools

# third party
import numpy as np

# first party
from src.operations import secrets
from src.utils.delphi_epidata import Epidata
from src.utils.epidate import EpiDate
from src.utils.epiweek import add_epiweeks, range_epiweeks
from src.utils.geo.locations import Locations


class DataSource(abc.ABC):
    """The interface by which all input data is provided."""

    @abc.abstractmethod
    def get_truth_locations(self):
        """Return a list of locations in which ground truth is available."""

    @abc.abstractmethod
    def get_sensor_locations(self):
        """Return a list of locations in which sensors are available."""

    @abc.abstractmethod
    def get_missing_locations(self, epiweek):
        """Return a tuple of locations which did not report on the given week."""

    @abc.abstractmethod
    def get_sensors(self):
        """Return a list of sensor names."""

    @abc.abstractmethod
    def get_weeks(self):
        """Return a list of weeks on which truth and sensors are both available."""

    @abc.abstractmethod
    def get_truth_value(self, epiweek, location):
        """Return ground truth (w)ILI."""

    @abc.abstractmethod
    def get_sensor_value(self, epiweek, location, name):
        """Return a sensor reading."""


class FluDataSource(DataSource):
    """The interface by which all input data is provided."""

    # the first epiweek for which we have ground truth ILI in all locations
    FIRST_DATA_EPIWEEK = 201040

    # all known sensors, past and present
    SENSORS = ['gft', 'ght', 'twtr', 'wiki', 'cdc', 'epic', 'sar3', 'arch']

    @staticmethod
    def new_instance():
        return FluDataSource(
            Epidata, FluDataSource.SENSORS, Locations.region_list)

    def __init__(self, epidata, sensors, locations):
        self.epidata = epidata
        self.sensors = sensors
        self.sensor_locations = locations

        # cache for prefetching bulk flu data
        self.cache = {}

    @functools.lru_cache(maxsize=1)
    def get_truth_locations(self):
        """Return a list of locations in which ground truth is available."""
        return Locations.region_list

    @functools.lru_cache(maxsize=1)
    def get_sensor_locations(self):
        """Return a list of locations in which sensors are available."""
        return self.sensor_locations

    @functools.lru_cache(maxsize=None)
    def get_missing_locations(self, epiweek):
        """Return a tuple of locations which did not report on the given week."""

        # only return missing atoms, i.e. locations that can't be further split
        atomic_locations = set(Locations.atom_list)

        available_locations = []
        for loc in atomic_locations:
            if self.get_truth_value(epiweek, loc) is None:
                # this atomic location didn't report (or it's a future week)
                continue
            available_locations.append(loc)

        if available_locations:
            return tuple(atomic_locations - set(available_locations))
        else:
            # no data is available, assume that all locations will be reporting
            return ()

    @functools.lru_cache(maxsize=1)
    def get_sensors(self):
        """Return a list of sensor names."""
        return self.sensors

    @functools.lru_cache(maxsize=1)
    def get_weeks(self):
        """Return a list of weeks on which truth and sensors are both available."""
        latest_week = self.get_most_recent_issue()
        week_range = range_epiweeks(
            FluDataSource.FIRST_DATA_EPIWEEK, latest_week, inclusive=True)
        return list(week_range)

    def get_truth_value(self, epiweek, location):
        """Return ground truth (w)ILI."""

        try:
            return self.cache['ilinet'][location][epiweek]
        except KeyError:
            # print('cache miss: get_truth_value', epiweek, location)
            auth = secrets.api.fluview
            response = self.epidata.fluview(location, epiweek, auth=auth)
            if response['result'] != 1:
                return self.add_to_cache('ilinet', location, epiweek, None)
            data = response['epidata'][0]
            if data['num_providers'] == 0:
                return self.add_to_cache('ilinet', location, epiweek, None)
            return self.add_to_cache('ilinet', location, epiweek, data['wili'])

    @functools.lru_cache(maxsize=None)
    def get_truth_values(self, epiweeks, location):
        """Return multiple ground truth (w)ILI."""
        try:
            return [self.cache['ilinet'][location][week] for week in epiweeks]
        except KeyError:
            auth = secrets.api.fluview
            response = self.epidata.fluview(location, epiweeks, auth=auth)
            if response['result'] != 1:
                return [self.add_to_cache('ilinet', location, week, None) for
                        week in epiweeks]

            valid_weeks = []
            for row in response['epidata']:
                valid_weeks.append(row['epiweek'])
            valid_weeks = np.array(valid_weeks)

            values = []
            for week in epiweeks:
                idx = list(np.where(valid_weeks == week))[0]
                if len(idx) > 0:
                    row = response['epidata'][idx[0]]
                    if row['num_providers'] == 0:
                        values.append(self.add_to_cache('ilinet', location,
                                                        row['epiweek'], None))
                        continue
                    values.append(
                        self.add_to_cache('ilinet', location, row['epiweek'],
                                          row['wili']))
                else:
                    values.append(
                        self.add_to_cache('ilinet', location, week, None))

        return values

    @functools.lru_cache(maxsize=None)
    def get_sensor_value(self, epiweek, location, name):
        """Return a sensor reading."""
        try:
            return self.cache[name][location][epiweek]
        except KeyError:
            # print('cache miss: get_sensor_value', epiweek, location, name)
            response = self.epidata.sensors(
                secrets.api.sensors, name, location, epiweek)
            if response['result'] != 1:
                return self.add_to_cache(name, location, epiweek, None)
            value = response['epidata'][0]['value']
            return self.add_to_cache(name, location, epiweek, value)

    @functools.lru_cache(maxsize=None)
    def get_sensor_values(self, epiweeks, location, name):
        """Return multiple sensor readings."""
        try:
            return [self.cache[name][location][week] for week in epiweeks]
        except KeyError:
            response = self.epidata.sensors(secrets.api.sensors, name, location,
                                            epiweeks)
            if response['result'] != 1:
                return [self.add_to_cache(name, location, week, None) for week
                        in epiweeks]
            valid_weeks = []
            for row in response['epidata']:
                valid_weeks.append(row['epiweek'])
            valid_weeks = np.array(valid_weeks)
            values = []
            for week in epiweeks:
                idx = list(np.where(valid_weeks == week))[0]
                if len(idx) > 0:
                    values.append(self.add_to_cache(name, location, week,
                                                    response['epidata'][idx[0]][
                                                        'value']))
                else:
                    values.append(self.add_to_cache(name, location, week, None))

        return values

    @functools.lru_cache(maxsize=1)
    def get_most_recent_issue(self):
        """Return the most recent epiweek for which FluView data is available."""
        ew2 = EpiDate.today().get_ew()
        ew1 = add_epiweeks(ew2, -9)
        response = self.epidata.fluview('nat', self.epidata.range(ew1, ew2))
        issues = [row['issue'] for row in self.epidata.check(response)]
        return max(issues)

    def multi_add_to_cache(self, name, location, epiweeks, values):
        """Add mulitple epiweek values to the cache."""
        if name not in self.cache:
            self.cache[name] = {}
        if location not in self.cache[name]:
            self.cache[name][location] = {}
        for week, value in zip(epiweeks, values):
            self.cache[name][location][week] = value

    def add_to_cache(self, name, location, epiweek, value):
        """Add the given value to the cache."""
        if name not in self.cache:
            self.cache[name] = {}
        if location not in self.cache[name]:
            self.cache[name][location] = {}
        self.cache[name][location][epiweek] = value
        return value
