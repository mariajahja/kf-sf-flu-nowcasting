"""
This file was created by the CMU DELPHI Epi-forecasting group.
https://github.com/cmu-delphi/

===============
=== Purpose ===
===============

Epiweek arithmetic and utilities.
"""


def split_epiweek(epiweek):
    """ return a (year, week) pair from this epiweek """
    return (epiweek // 100, epiweek % 100)


def join_epiweek(year, week):
    """ return an epiweek from the (year, week) pair """
    return year * 100 + week


def check_epiweek(*epiweeks, fatal=True):
    """ return True if the epiweek is valid, otherwise raise Exception """
    for epiweek in epiweeks:
        year, week = split_epiweek(epiweek)
        if not 1 <= week <= get_num_weeks(year):
            if fatal:
                raise Exception('invalid epiweek: epiweek=%d' % epiweek)
            else:
                return False
    return True


def get_num_weeks(year):
    """ return the number of epiweeks in the year """
    if not 1900 <= year < 2100:
        raise Exception('not sure how many epiweeks: year=%d' % year)
    elif (year % 28) in [4, 9, 15, 20, 26]:
        return 53
    else:
        return 52


def add_epiweeks(epiweek, i):
    """ return the epiweek plus (or minus) the number of weeks """
    check_epiweek(epiweek)
    year, week = split_epiweek(epiweek)
    while i > 0:
        weeks_in_year = get_num_weeks(year)
        dw = min(i, weeks_in_year - week)
        i -= dw
        week += dw
        if i > 0 and week == weeks_in_year:
            i -= 1
            year += 1
            week = 1
    while i < 0:
        dw = min(week - 1, -i)
        i += dw
        week -= dw
        if i < 0 and week == 1:
            i += 1
            year -= 1
            week = get_num_weeks(year)
    return join_epiweek(year, week)


def get_season(epiweek, offseason=lambda x: (None, None)):
    """ return the epiweek range of the flu season containing the epiweek """
    check_epiweek(epiweek)
    year, week = split_epiweek(epiweek)
    if week <= 20:
        return (join_epiweek(year - 1, 40), join_epiweek(year, 20))
    elif week >= 40:
        return (join_epiweek(year, 40), join_epiweek(year + 1, 20))
    else:
        return offseason(epiweek)


def delta_epiweeks(ew1, ew2):
    """ return the number of weeks between the two epiweeks """
    check_epiweek(ew1, ew2)
    y1, w1 = split_epiweek(ew1)
    y2, w2 = split_epiweek(ew2)
    num = 0
    while y2 > y1:
        num += get_num_weeks(y2 - 1)
        y2 -= 1
    while y2 < y1:
        num -= get_num_weeks(y2)
        y2 += 1
    num += w2 - w1
    return num


def range_epiweeks(start, stop=None, inclusive=False, num=None):
    """
    an epiweek generator function
    exactly one of "stop" and "num" must be specified
    like Python's built-in "(x)range" function, "stop" is exclusive unless
      otherwise specified with the "inclusive" parameter
    """
    if (stop is None) == (num is None):
        raise Exception('Exactly one of "stop" and "num" must be specified')
    if num is None:
        num = delta_epiweeks(start, stop)
        if inclusive:
            if stop >= start:
                num += 1
            else:
                num -= 1
    while num > 0:
        yield start
        num -= 1
        start = add_epiweeks(start, 1)
    while num < 0:
        yield start
        num += 1
        start = add_epiweeks(start, -1)
