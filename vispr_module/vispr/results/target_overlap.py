# coding: utf-8
from __future__ import absolute_import, division, print_function

__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster, Liu lab"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


import json
from itertools import combinations
from operator import itemgetter


def overlap(*targets):
    isect = set(targets[0])
    for other in targets[1:]:
        isect &= other
    return isect


def overlaps(order, **targets):
    """
    Arguments
    order   -- 1: single condition, 2: overlap of 3 conditions, 3: overlap of 3 conditions...
    targets -- labels and targets to compare
    """
    for c in combinations(targets.items(), order):
        isect = overlap(*map(itemgetter(1), c))
        labels = list(map(itemgetter(0), c))
        yield labels, len(isect)


def plot_overlap_venn(**targets):
    data = []
    for s in range(1, len(targets) + 1):
        for labels, isect in overlaps(s, **targets):
            data.append({"sets": labels, "size": isect})
    return json.dumps(data)
