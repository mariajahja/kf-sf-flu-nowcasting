"""
This file was created by the CMU DELPHI Epi-forecasting group.
https://github.com/cmu-delphi/

===============
=== Purpose ===
===============

Date arithmetic that's epiweek-aware. Based on EpiVis's Date class at:
https://github.com/cmu-delphi/www-epivis/blob/master/site/js/epivis.js
"""

# standard library
import datetime

# first party
from src.utils.epiweek import get_num_weeks


class EpiDate:
    DAYS_PER_MONTH = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    CUMULATIVE_DAYS_PER_MONTH = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273,
                                 304, 334]
    DAY_OF_WEEK_TABLE = [0, 3, 2, 5, 0, 3, 5, 1, 4, 6, 2, 4]
    DAY_NAMES_LONG = ['Sunday', 'Monday', 'Tuesday', 'Wednesday', 'Thursday',
                      'Friday', 'Saturday']
    DAY_NAMES_SHORT = ['Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat']
    MONTH_NAMES_LONG = ['January', 'February', 'March', 'April', 'May', 'June',
                        'July', 'August', 'September', 'October', 'November',
                        'December']
    MONTH_NAMES_SHORT = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug',
                         'Sep', 'Oct', 'Nov', 'Dec']

    def __init__(self, year, month, day):
        if year < 1 or month < 1 or month > 12 or day < 1 or day > \
                EpiDate.DAYS_PER_MONTH[month - 1] + (
                1 if month == 2 and EpiDate._is_leap_year(year) else 0):
            raise Exception('invalid date: %d/%d/%d' % (year, month, day))
        self.year = year
        self.month = month
        self.day = day

    def is_leap_year(self):
        return EpiDate._is_leap_year(self.year)

    def get_index(self):
        return EpiDate._get_index(self.year, self.month, self.day)

    def get_day_of_week(self):
        return EpiDate._get_day_of_week(self.year, self.month, self.day)

    def get_year(self):
        return self.year

    def get_month(self):
        return self.month

    def get_day(self):
        return self.day

    def get_ew_year(self):
        this_date = self.get_index()
        first_date1 = EpiDate._get_index(self.year, 1,
                                         4) - EpiDate._get_day_of_week(
            self.year, 1, 4)
        if this_date < first_date1:
            y = -1
        else:
            first_date2 = EpiDate._get_index(self.year + 1, 1,
                                             4) - EpiDate._get_day_of_week(
                self.year + 1, 1, 4)
            y = 0 if this_date < first_date2 else 1
        return self.year + y

    def get_ew_week(self):
        this_date = self.get_index()
        first_date1 = EpiDate._get_index(self.year, 1,
                                         4) - EpiDate._get_day_of_week(
            self.year, 1, 4)
        if this_date < first_date1:
            first_date = EpiDate._get_index(self.year - 1, 1,
                                            4) - EpiDate._get_day_of_week(
                self.year - 1, 1, 4)
        else:
            first_date2 = EpiDate._get_index(self.year + 1, 1,
                                             4) - EpiDate._get_day_of_week(
                self.year + 1, 1, 4)
            first_date = first_date1 if this_date < first_date2 else first_date2
        return ((this_date - first_date) // 7) + 1

    def get_ew(self):
        return self.get_ew_year() * 100 + self.get_ew_week()

    def add_days(self, num):
        return EpiDate.from_index(self.get_index() + num)

    def add_weeks(self, num):
        return self.add_days(num * 7)

    def add_months(self, num):
        m = self.year * 12 + (self.month - 1) + num
        year = m // 12
        month = (m % 12) + 1
        leap_day = 1 if month == 2 and EpiDate._is_leap_year(year) else 0
        max_day = EpiDate.DAYS_PER_MONTH[month - 1] + leap_day
        return EpiDate(year, month, min(self.day, max_day))

    def add_years(self, num):
        return self.add_months(num * 12)

    def __str__(self):
        return '%04d-%02d-%02d' % (self.year, self.month, self.day)

    @staticmethod
    def get_day_name(day, short=False):
        names = EpiDate.DAY_NAMES_SHORT if short else EpiDate.DAY_NAMES_LONG
        return names[day]

    @staticmethod
    def get_month_name(month, short=False):
        names = EpiDate.MONTH_NAMES_SHORT if short else EpiDate.MONTH_NAMES_LONG
        return names[month - 1]

    @staticmethod
    def today():
        today = datetime.datetime.now()
        return EpiDate(today.year, today.month, today.day)

    @staticmethod
    def from_string(str):
        if len(str) == 8:
            return EpiDate(int(str[0:4]), int(str[4:6]), int(str[6:8]))
        elif len(str) == 10:
            return EpiDate(int(str[0:4]), int(str[5:7]), int(str[8:10]))
        else:
            raise Exception('expected YYYYMMDD or YYYY_MM_DD')

    @staticmethod
    def from_index(index):
        x = index
        year, index = index // 365, index % 365
        leaps = (year // 4) - (year // 100) + (year // 400)
        index = index - leaps + year * 365
        year = (index // 365) + 1
        y = (year + 399) % 400
        year_offset = (((year - 1) // 400) * 146097) + (y * 365) + (y // 4) - (
                y // 100) + (y // 400)
        if x - year_offset == 365 and EpiDate._is_leap_year(year + 1):
            year += 1
            year_offset += 365
        index = x - year_offset
        m = 1
        leap = 1 if EpiDate._is_leap_year(year) else 0
        while m < 12 and EpiDate.CUMULATIVE_DAYS_PER_MONTH[m] + (
                leap if m >= 2 else 0) <= index:
            m += 1
        index -= EpiDate.CUMULATIVE_DAYS_PER_MONTH[m - 1] + (
            1 if m > 2 and EpiDate._is_leap_year(year) else 0)
        return EpiDate(year, m, index + 1)

    @staticmethod
    def from_epiweek(year, week):
        if year < 1 or week < 1 or week > get_num_weeks(year):
            raise Exception('invalid year or week')
        date = EpiDate(year, 7, 1)
        while date.get_ew_week() < week:
            date = date.add_weeks(+1)
        while date.get_ew_week() > week:
            date = date.add_weeks(-1)
        while date.get_day_of_week() < 3:
            date = date.add_days(+1)
        while date.get_day_of_week() > 3:
            date = date.add_days(-1)
        return date

    @staticmethod
    def _is_leap_year(year):
        return (year % 4 == 0 and year % 100 != 0) or year % 400 == 0

    @staticmethod
    def _get_index(year, month, day):
        y = (year + 399) % 400
        year_offset = (((year - 1) // 400) * 146097) + (y * 365) + (y // 4) - (
                y // 100) + (y // 400)
        month_offset = EpiDate.CUMULATIVE_DAYS_PER_MONTH[month - 1] + (
            1 if month > 2 and EpiDate._is_leap_year(year) else 0)
        day_offset = day - 1
        return year_offset + month_offset + day_offset

    @staticmethod
    def _get_day_of_week(year, month, day):
        y = year - (1 if month < 3 else 0)
        return (y + (y // 4) - (y // 100) + (y // 400) +
                EpiDate.DAY_OF_WEEK_TABLE[month - 1] + day) % 7
