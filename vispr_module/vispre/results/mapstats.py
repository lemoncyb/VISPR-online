# coding: utf-8
from __future__ import absolute_import, division, print_function

__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster, Liu lab"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import os

import numpy as np
import pandas as pd
from flask import render_template

from ..results.common import AbstractResults, templates

class Results(AbstractResults):
    def __init__(self, dataframe):
        super(Results, self).__init__(dataframe)
        self.df.columns = self.df.columns.str.lower()
        self.df["file"] = self.df["file"].apply(os.path.basename).apply(lambda filename: os.path.splitext(os.path.splitext(filename)[0])[0])
        self.df.sort_values("file", inplace=True)
        self.has_replicates = self.df["label"].value_counts()[0] > 1

    def plot_mapstats(self):
        data = self.df.loc[:, ["file", "label", "reads", "mapped"]]
        data["unmapped_percentage"] = ((self.df["reads"] - self.df["mapped"]) / self.df["reads"]).apply("{:.1%}".format)

        width = 20 * data.shape[0]
        plt = templates.get_template("plots/mapstats.json").render(
                               data=data.to_json(orient="records"),
                               width=width,
                               padding=0.1 if self.has_replicates else 0)
        return plt

    def plot_zerocounts(self):
        data = self.df.loc[:, ["file", "label"]]
        data["zerocounts"] = np.log10(self.df["zerocounts"])

        width = 20 * data.shape[0]
        return templates.get_template("plots/zerocounts.json").render(
                               data=data.to_json(orient="records"),
                               width=width,
                               padding=0.1 if self.has_replicates else 0)

    def plot_gini_index(self):
        data = self.df[["file", "label", "giniindex"]]

        width = 20 * data.shape[0]
        return templates.get_template("plots/gini_index.json").render(
                               data=data.to_json(orient="records"),
                               width=width,
                               padding=0.1 if self.has_replicates else 0)
