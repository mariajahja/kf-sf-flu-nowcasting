"""
This file was created by the CMU DELPHI Epi-forecasting group.
https://github.com/cmu-delphi/

===============
=== Purpose ===
===============

Encodes the hierarchy of US political divisions.

This file, together with populations.py, replaces the static data portion of
state_info.py.

The location names used in this file match FluView names as specified in
fluview_locations.py of delphi-epidata.


===================
=== Explanation ===
===================

Although intended to be a more or less general-purpose description of the
various US geopolitical divisions, for all practical purposes the data in this
file corresponds to the FluView perspective of the world.

In this perspective, the US is a hierarchy where regions at any given level are
composed of smaller regions at a lower level. Notably, it may be possible to
subdivide a given region into multiple distinct sets of smaller regions.
However, the set of locations in any given subdivision fully covers and spans
the region being subdivided. In other words, there are never any gaps.

The root of the hierarchy is the national region (shortened to "nat") which
represents the entire US, including many of its territories. Each lower layer
of the hierarchy consists of smaller regions which combine together to form the
national region.

The leaves of the hierarchy are called "atoms" and have no further subdivisions
-- at least, not from a FluView perspective. These are typically US states,
although they also include some state fragments, territories, and cities.

By convention, the the middle layers of the hierarchy are collectively called
"regions". This includes, for example, the ten HHS regions as one subdivision
of national and the nine Census divisions as another. Each of the HHS and
Census regions is in turn made up of atoms -- mostly states, with a few
exceptions.
"""


class Locations:
    """Encodes the hierarchy of US political divisions."""

    # atomic regions for FluView data (regions containing only themselves)
    atom_list = [
        # entire states
        'ak', 'al', 'ar', 'az', 'ca', 'co', 'ct', 'de', 'fl', 'ga', 'hi', 'ia',
        'id', 'il', 'in', 'ks', 'ky', 'la', 'ma', 'md', 'me', 'mi', 'mn', 'mo',
        'ms', 'mt', 'nc', 'nd', 'ne', 'nh', 'nj', 'nm', 'nv', 'oh', 'ok', 'or',
        'pa', 'ri', 'sc', 'sd', 'tn', 'tx', 'ut', 'va', 'vt', 'wa', 'wi', 'wv',
        'wy',
        # state fragments
        'ny_minus_jfk',
        # territories
        'dc', 'pr', 'vi',
        # cities
        'jfk',
    ]
    atom_map = {a: [a] for a in atom_list}

    # national, HHS, and Census regions in terms of atoms
    nat_list = ['nat']
    nat_map = dict(zip(nat_list, [atom_list]))

    hhs_list = ['hhs%d' % i for i in range(1, 11)]
    hhs_map = dict(zip(hhs_list, [
        ['ct', 'ma', 'me', 'nh', 'ri', 'vt'],
        ['jfk', 'nj', 'ny_minus_jfk', 'pr', 'vi'],
        ['dc', 'de', 'md', 'pa', 'va', 'wv'],
        ['al', 'fl', 'ga', 'ky', 'ms', 'nc', 'sc', 'tn'],
        ['il', 'in', 'mi', 'mn', 'oh', 'wi'],
        ['ar', 'la', 'nm', 'ok', 'tx'],
        ['ia', 'ks', 'mo', 'ne'],
        ['co', 'mt', 'nd', 'sd', 'ut', 'wy'],
        ['az', 'ca', 'hi', 'nv'],
        ['ak', 'id', 'or', 'wa'],
    ]))

    cen_list = ['cen%d' % i for i in range(1, 10)]
    cen_map = dict(zip(cen_list, [
        ['ct', 'ma', 'me', 'nh', 'ri', 'vt'],
        ['jfk', 'nj', 'ny_minus_jfk', 'pa', 'pr', 'vi'],
        ['il', 'in', 'mi', 'oh', 'wi'],
        ['ia', 'ks', 'mn', 'mo', 'nd', 'ne', 'sd'],
        ['dc', 'de', 'fl', 'ga', 'md', 'nc', 'sc', 'va', 'wv'],
        ['al', 'ky', 'ms', 'tn'],
        ['ar', 'la', 'ok', 'tx'],
        ['az', 'co', 'id', 'mt', 'nm', 'nv', 'ut', 'wy'],
        ['ak', 'ca', 'hi', 'or', 'wa'],
    ]))

    # New York state combines the "ny_minus_jfk" fragment with the "jfk" city
    ny_state_list = ['ny']
    ny_state_map = {ny_state_list[0]: ['jfk', 'ny_minus_jfk']}

    # collections of all known locations
    region_list = nat_list + hhs_list + cen_list + ny_state_list + atom_list
    region_map = {}
    region_map.update(nat_map)
    region_map.update(hhs_map)
    region_map.update(cen_map)
    region_map.update(ny_state_map)
    region_map.update(atom_map)
