"""
This file was created by the CMU DELPHI Epi-forecasting group.
https://github.com/cmu-delphi/

===============
=== Purpose ===
===============

Contains imputed state-level weights for past seasons. This is helpful because
although CDC publishes detailed ILINet data, they don't publish the population
values used in computing wILI. Rather than using population estimates, we use
the (relative) population weights imputed from the data.

This file should be updated annually on or after data has been released for
epiweek 44 (i.e. Friday of epiweek 45, whicih is typically the second week of
November).


===================
=== Methodology ===
===================

The weights in this file were imputed using some manual reverse engineering
followed by sensor fusion (or regression, depending on how you look at it).
There are two important points to highlight:

  1. It appears to be the case that CDC uses a fixed set of population numbers
  within each season. These values are updated each year on epiweek 40 and
  remain constant through epiweek 39 of the following year.

  2. Occasionally a location will report an ILI of zero, which implies that
  num_ili is zero. But this doesn't imply that num_providers is also zero. It
  appears to the the case that ILI of zero is included in wILI iff
  num_providers is nonzero. This means that the composition of each region is
  subject to change each epiweek, depending on whether there are any reporting
  providers in each location on that week.

Given the above observations, the goal is to determine the (relative) weights
used to compute national wILI as a function of state-level ILI, separately for
each season (starting in 2010, since that's when state-level data publicly
begins).

Here's a first approach. Define column vector Y as national wILI on each of the
52 (or 53) weeks of the season. Define X to be the matrix with columns
representing per-state (or territory or city) ILI on the same weeks. Solve for
weight vector B such that X * B = Y. This can be done in the spirit of
regression as: B = (X^T * X)^{-1} * X^T * Y. As a sanity check, the entries of
B should sum to one.

However, there are two problems with the above formulation. First, there are
sporadic missing values, caused by zero reporters (not the same as zero ILI) in
some locations on some weeks. Second, the system is underdetermined because the
number of weeks is 52 (or 53), but the number of locations is currently 54 and
likely to increase in the future.

We solve the first by building X using regional wILI in place of state ILI in
locations and on weeks where there are no providers. We solve the second by
separating the imputation into two stages. We first find the weights of the HHS
regions relative to national, and then we find the weights of states relative
to their corresponding HHS regions. Finally, we combine the state-to-regional
and regional-to-national weights to get the weights relating state ILI to
national wILI -- the data in this file.

Since the weights must necessarily sum to one, we compute these weights via
sensor fusion out of convenience. Define H as a column vector of ones, Z as the
matrix where Y is subtracted from each column of X (i.e. Z = X - Y * H^T), and
R as the empirical covariance matrix of Z (i.e. R ∝ Z^T * Z). The weight vector
B is computed as: B = R^{-1} * H * (H^T * R^{-1} * H)^{-1}.

Since the size of each group is smaller (e.g. ~5 states per region instead of
all 54 state-level locations), fewer weeks of data are needed. The largest
group is the 10 HHS regions (no single region has that many states), and so 10
weeks of data are needed to exactly impute all weights. (This could actually be
improved by using the 9 Census divisions instead of the 10 HHS regions.)

By exploiting backfill, it's possible to learn the weights for all state-level
locations just four weeks into the season, i.e. upon publication of issue w43.
(At that point, there will be 4 versions of w40, 3 of w41, 2 of w42 and 1 of
w43 -- 10 in total.)


===============
=== Results ===
===============

The imputed weights were used to compute wILI for all regions and on all weeks
from 2010w40 through 2018w01. Computed and reported wILI differ by less than
0.0015 ILI, with a median difference on the order of 0.00001. For this reason,
and in the interest of brevity, values in this file are conservatively rounded
to eight decimal places. For reference, reported wILI is rounded to between 3
and 5 decimal places, depending on data source and reporting date.

A cursory investigation suggests that these imputed weights differ from the
fixed population estimates we previously used by, in some cases, up to 10%.
This has potential to meaningfully improve wILI nowcasts, especially in smaller
regions.
"""

population_weights = {
    2010: {
        'ak': 0.00227518, 'al': 0.01533745, 'ar': 0.00941176, 'az': 0.02148400,
        'ca': 0.12039274, 'co': 0.01636543, 'ct': 0.01146016, 'dc': 0.00195314,
        'de': 0.00288299, 'fl': 0.06038341, 'ga': 0.03201648, 'hi': 0.00421877,
        'ia': 0.00979696, 'id': 0.00503527, 'il': 0.04205408, 'in': 0.02092240,
        'jfk': 0.02733469, 'ks': 0.00918118, 'ky': 0.01405184, 'la': 0.01463192,
        'ma': 0.02147752, 'md': 0.01856371, 'me': 0.00429419, 'mi': 0.03247482,
        'mn': 0.01715415, 'mo': 0.01950239, 'ms': 0.00961546, 'mt': 0.00317556,
        'nc': 0.03055618, 'nd': 0.00210670, 'ne': 0.00585185, 'nh': 0.00431464,
        'nj': 0.02836354, 'nm': 0.00654612, 'nv': 0.00860921,
        'ny_minus_jfk': 0.03631738, 'oh': 0.03759897, 'ok': 0.01200984,
        'or': 0.01246166, 'pa': 0.04105459, 'ri': 0.00343073, 'sc': 0.01485699,
        'sd': 0.00264588, 'tn': 0.02050899, 'tx': 0.08072283, 'ut': 0.00906917,
        'va': 0.02567430, 'vt': 0.00202532, 'wa': 0.02170791, 'wi': 0.01841977,
        'wv': 0.00592713, 'wy': 0.00177260,
    },
    2011: {
        'ak': 0.00229505, 'al': 0.01547462, 'ar': 0.00943623, 'az': 0.02070624,
        'ca': 0.12068005, 'co': 0.01632952, 'ct': 0.01153898, 'dc': 0.00195012,
        'de': 0.00291042, 'fl': 0.06087038, 'ga': 0.03136422, 'hi': 0.00440658,
        'ia': 0.00987440, 'id': 0.00506529, 'il': 0.04159759, 'in': 0.02102140,
        'jfk': 0.02639945, 'ks': 0.00924801, 'ky': 0.01404776, 'la': 0.01467032,
        'ma': 0.02113917, 'md': 0.01871206, 'me': 0.00428860, 'mi': 0.03204299,
        'mn': 0.01719550, 'mo': 0.01941239, 'ms': 0.00960674, 'mt': 0.00321277,
        'nc': 0.03087160, 'nd': 0.00218381, 'ne': 0.00591985, 'nh': 0.00425016,
        'nj': 0.02839098, 'nm': 0.00666356, 'nv': 0.00874816,
        'ny_minus_jfk': 0.03617721, 'oh': 0.03740244, 'ok': 0.01213960,
        'or': 0.01237936, 'pa': 0.04116837, 'ri': 0.00339827, 'sc': 0.01497568,
        'sd': 0.00264368, 'tn': 0.02054620, 'tx': 0.08137335, 'ut': 0.00897433,
        'va': 0.02593179, 'vi': 0.00034350, 'vt': 0.00202018, 'wa': 0.02172891,
        'wi': 0.01843688, 'wv': 0.00600521, 'wy': 0.00183008,
    },
    2012: {
        'ak': 0.00231819, 'al': 0.01540851, 'ar': 0.00942543, 'az': 0.02079732,
        'ca': 0.12092341, 'co': 0.01641994, 'ct': 0.01148834, 'dc': 0.00198275,
        'de': 0.00291038, 'fl': 0.06114108, 'ga': 0.03148930, 'hi': 0.00441076,
        'ia': 0.00982549, 'id': 0.00508404, 'il': 0.04128530, 'in': 0.02090667,
        'jfk': 0.02644989, 'ks': 0.00921252, 'ky': 0.01401796, 'la': 0.01467680,
        'ma': 0.02113535, 'md': 0.01869878, 'me': 0.00426125, 'mi': 0.03168319,
        'mn': 0.01714633, 'mo': 0.01928559, 'ms': 0.00955571, 'mt': 0.00320334,
        'nc': 0.03098006, 'nd': 0.00219478, 'ne': 0.00591215, 'nh': 0.00422915,
        'nj': 0.02829857, 'nm': 0.00668008, 'nv': 0.00873694,
        'ny_minus_jfk': 0.03599504, 'oh': 0.03703687, 'ok': 0.01216368,
        'or': 0.01241930, 'pa': 0.04088358, 'ri': 0.00337286, 'sc': 0.01501189,
        'sd': 0.00264458, 'tn': 0.02054316, 'tx': 0.08236755, 'ut': 0.00904064,
        'va': 0.02597622, 'vi': 0.00035180, 'vt': 0.00200973, 'wa': 0.02190810,
        'wi': 0.01832391, 'wv': 0.00595245, 'wy': 0.00182329,
    },
    2013: {
        'ak': 0.00230241, 'al': 0.01517830, 'ar': 0.00928298, 'az': 0.02062798,
        'ca': 0.11974465, 'co': 0.01633128, 'ct': 0.01130065, 'dc': 0.00199037,
        'de': 0.00288681, 'fl': 0.06080683, 'ga': 0.03122503, 'hi': 0.00438272,
        'ia': 0.00967646, 'id': 0.00502290, 'il': 0.04052832, 'in': 0.02057813,
        'jfk': 0.02624223, 'ks': 0.00908388, 'ky': 0.01378813, 'la': 0.01448562,
        'ma': 0.02091874, 'md': 0.01852285, 'me': 0.00418369, 'mi': 0.03111046,
        'mn': 0.01693229, 'mo': 0.01895510, 'ms': 0.00939586, 'mt': 0.00316437,
        'nc': 0.03069718, 'nd': 0.00220253, 'ne': 0.00584064, 'nh': 0.00415701,
        'nj': 0.02790382, 'nm': 0.00656480, 'nv': 0.00868439,
        'ny_minus_jfk': 0.03536077, 'oh': 0.03633856, 'ok': 0.01200804,
        'or': 0.01227398, 'pa': 0.04017603, 'pr': 0.01154317, 'ri': 0.00330567,
        'sc': 0.01486956, 'sd': 0.00262349, 'tn': 0.02032259, 'tx': 0.08202735,
        'ut': 0.00898886, 'va': 0.02576680, 'vi': 0.00033496, 'vt': 0.00197035,
        'wa': 0.02170973, 'wi': 0.01802582, 'wv': 0.00584026, 'wy': 0.00181457,
    },
    2014: {
        'ak': 0.00229826, 'al': 0.01511231, 'ar': 0.00925232, 'az': 0.02071714,
        'ca': 0.11984244, 'co': 0.01647139, 'ct': 0.01124367, 'dc': 0.00202117,
        'de': 0.00289444, 'fl': 0.06113065, 'ga': 0.03123955, 'hi': 0.00438962,
        'ia': 0.00966159, 'id': 0.00504010, 'il': 0.04027598, 'in': 0.02054438,
        'jfk': 0.02628095, 'ks': 0.00904742, 'ky': 0.01374139, 'la': 0.01446115,
        'ma': 0.02092603, 'md': 0.01853664, 'me': 0.00415316, 'mi': 0.03093930,
        'mn': 0.01694700, 'mo': 0.01889603, 'ms': 0.00935160, 'mt': 0.00317376,
        'nc': 0.03078912, 'nd': 0.00226162, 'ne': 0.00584164, 'nh': 0.00413793,
        'nj': 0.02782429, 'nm': 0.00651961, 'nv': 0.00872335,
        'ny_minus_jfk': 0.03515874, 'oh': 0.03617563, 'ok': 0.01203859,
        'or': 0.01228678, 'pa': 0.03993758, 'pr': 0.01130270, 'ri': 0.00328770,
        'sc': 0.01492810, 'sd': 0.00264159, 'tn': 0.02030914, 'tx': 0.08268967,
        'ut': 0.00906940, 'va': 0.02582641, 'vi': 0.00033273, 'vt': 0.00195921,
        'wa': 0.02179514, 'wi': 0.01795467, 'wv': 0.00579754, 'wy': 0.00182168,
    },
    2015: {
        'ak': 0.00222549, 'al': 0.01465046, 'ar': 0.00896152, 'az': 0.02033757,
        'ca': 0.11723059, 'co': 0.01617926, 'ct': 0.01086707, 'dc': 0.00199051,
        'de': 0.00282642, 'fl': 0.06009969, 'ga': 0.03050535, 'hi': 0.00428888,
        'ia': 0.00938773, 'id': 0.00493738, 'il': 0.03891432, 'in': 0.01993020,
        'jfk': 0.02565219, 'ks': 0.00877424, 'ky': 0.01333317, 'la': 0.01404697,
        'ma': 0.02038083, 'md': 0.01805455, 'me': 0.00401881, 'mi': 0.02993998,
        'mn': 0.01648714, 'mo': 0.01832052, 'ms': 0.00904533, 'mt': 0.00309194,
        'nc': 0.03004141, 'nd': 0.00223386, 'ne': 0.00568488, 'nh': 0.00400892,
        'nj': 0.02700323, 'nm': 0.00630083, 'nv': 0.00857763,
        'ny_minus_jfk': 0.05965562, 'oh': 0.03502746, 'ok': 0.01171597,
        'or': 0.01199323, 'pa': 0.03862992, 'pr': 0.01072014, 'ri': 0.00318804,
        'sc': 0.01459964, 'sd': 0.00257716, 'tn': 0.01978616, 'tx': 0.08143924,
        'ut': 0.00888997, 'va': 0.02515309, 'vi': 0.00032150, 'vt': 0.00189317,
        'wa': 0.02133134, 'wi': 0.01739496, 'wv': 0.00558998, 'wy': 0.00176455,
    },
    2016: {
        'ak': 0.00227206, 'al': 0.01495049, 'ar': 0.00916311, 'az': 0.02100899,
        'ca': 0.12044360, 'co': 0.01679052, 'ct': 0.01104830, 'dc': 0.00206836,
        'de': 0.00291047, 'fl': 0.06237330, 'ga': 0.03143072, 'hi': 0.00440484,
        'ia': 0.00961212, 'id': 0.00509207, 'il': 0.03957187, 'in': 0.02036971,
        'jfk': 0.02630960, 'ks': 0.00895905, 'ky': 0.01361560, 'la': 0.01437120,
        'ma': 0.02090476, 'md': 0.01848052, 'me': 0.00409015, 'mi': 0.03053260,
        'mn': 0.01689232, 'mo': 0.01871922, 'ms': 0.00920701, 'mt': 0.00317852,
        'nc': 0.03090132, 'nd': 0.00232914, 'ne': 0.00583455, 'nh': 0.00409392,
        'nj': 0.02756388, 'nm': 0.00641636, 'nv': 0.00889466,
        'ny_minus_jfk': 0.03460203, 'oh': 0.03573599, 'ok': 0.01203489,
        'or': 0.01239663, 'pa': 0.03939073, 'pr': 0.01069004, 'ri': 0.00324996,
        'sc': 0.01506526, 'sd': 0.00264163, 'tn': 0.02030874, 'tx': 0.08452051,
        'ut': 0.00921879, 'va': 0.02579260, 'vi': 0.00032227, 'vt': 0.00192612,
        'wa': 0.02206222, 'wi': 0.01775927, 'wv': 0.00567395, 'wy': 0.00180349,
    },
    2017: {
        'ak': 0.00226985, 'al': 0.01488923, 'ar': 0.00914813, 'az': 0.02121906,
        'ca': 0.12016114, 'co': 0.01697014, 'ct': 0.01094766, 'dc': 0.00208648,
        'de': 0.00291668, 'fl': 0.06310502, 'ga': 0.03156409, 'hi': 0.00437360,
        'ia': 0.00959975, 'id': 0.00514951, 'il': 0.03919113, 'in': 0.02030593,
        'jfk': 0.02613827, 'ks': 0.00890329, 'ky': 0.01358380, 'la': 0.01433368,
        'ma': 0.02085150, 'md': 0.01842029, 'me': 0.00407546, 'mi': 0.03039408,
        'mn': 0.01689905, 'mo': 0.01865930, 'ms': 0.00914998, 'mt': 0.00319311,
        'nc': 0.03106493, 'nd': 0.00232165, 'ne': 0.00584048, 'nh': 0.00408617,
        'nj': 0.02738362, 'nm': 0.00637229, 'nv': 0.00900122,
        'ny_minus_jfk': 0.03431267, 'oh': 0.03555613, 'ok': 0.01200884,
        'or': 0.01252378, 'pa': 0.03913594, 'pr': 0.01044379, 'ri': 0.00323374,
        'sc': 0.01518871, 'sd': 0.00265086, 'tn': 0.02036199, 'tx': 0.08529818,
        'ut': 0.00934570, 'va': 0.02575484, 'vi': 0.00031504, 'vt': 0.00191183,
        'wa': 0.02229756, 'wi': 0.01769118, 'wv': 0.00560636, 'wy': 0.00179330,
    },
}

first_season, last_season = min(population_weights), max(population_weights)


def get_population_weight(location, season=last_season):
    """
    Return the population weight of the given location, relative to the US
    nationally.

    inputs:
      location: A FluView location, which is case-sensitive and is expected to be
        lower case. The locations available generally represent the finest
        geographic resolution available for FluView data; namely, states,
        territories, and cities (aka "atoms"). Regional weights are not
        precomputed here because they differ by epiweek, depending on which
        locations are reporting (i.e. whether num_providers is nonzero).
      season (optional): The year containing epiweek 40. For example, the 2017
        season spans 2017w40--2018w39. By default, the most recent data is used.

    output:
      - the fraction of the US population contained within the given location
    """
    season = max(min(season, last_season), first_season)
    return population_weights[season][location]


def get_population(location, season=last_season):
    """
    Return the approximate population of the given location. The returned value
    is rounded to the nearest integer and is based on the assumption that the US
    population is a fixed 325 million.

    inputs: see `get_population_weight`

    output:
      - the approimate population of the given location
    """
    return round(get_population_weight(location, season) * 3.25e8)
