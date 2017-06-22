# coding: utf-8
from __future__ import absolute_import, division, print_function

__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster, Liu lab"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


import pandas as pd
from jinja2 import Environment, FileSystemLoader
from config import TEMPLATES_DATA_DIR, COUNT_NORMALIZED_HEADER


templates = Environment(loader=FileSystemLoader(TEMPLATES_DATA_DIR))


try:
    from functools import lru_cache
except ImportError:
    # dummy cache, i.e. Python 2 version will be a bit slower.
    def lru_cache():
        def dummy(func):
            return func

        return dummy


class AbstractResults(object):
    def __init__(self, dataframe, samples):
        """
        Arguments:
        dataframe -- a pandas data frame or its path consisting of per gene results as produced by MAGeCK
        """
        cols = samples + COUNT_NORMALIZED_HEADER # construct selected columns in count_normalized file
        self.df = pd.read_table(dataframe, usecols=cols, na_filter=False, low_memory=False)

    def __getitem__(self, slice):
        return self.df.__getitem__(slice)
