# coding: utf-8
from __future__ import absolute_import, division, print_function

__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster, Liu lab"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


from flask import render_template
import pandas as pd
import numpy as np

from ..results.common import lru_cache, AbstractResults, templates


class Results(AbstractResults):
    """Keep and display target results."""

    def __init__(self, dataframe, table_filter=None):
        """
        Arguments

        dataframe -- path to file containing MAGeCK target (gene) summary. Alternatively, a dataframe.
        table_filter  -- an optional function (taking a results row and returning bool) to filter the stored results. This will not affect p-value distribution plots. Default: None
        """
        self.df = dataframe
        if self.df.index.duplicated().any():
            raise ValueError("Target results contain duplicated gene names.")
        self.df.sort_values("p-value", inplace=True)
        self.df.reset_index(drop=True, inplace=True)
        self.df["log10-p-value"] = -np.log10(self.df["p-value"])
        self.df["idx"] = self.df.index
        self.df.index = self.df["target"]

        pvals = self.df["log10-p-value"]
        noninf = pvals.replace([np.inf, -np.inf], np.nan).dropna()
        if noninf.empty:
            vmax = 0
            vmin = 0
        else:
            vmax, vmin = noninf.max(), noninf.min()
        pval_cdf = pvals.replace(np.inf, vmax) \
                        .replace(-np.inf, vmin)
        pval_cdf = pval_cdf.value_counts(normalize=True, sort=False, bins=1000).cumsum()
        pval_cdf.index = np.maximum(0, pval_cdf.index)
        self.pval_cdf = pd.DataFrame({"p-value": pval_cdf.index, "cdf": pval_cdf})

        edges = np.arange(0, 1.01, 0.05)
        counts, _ = np.histogram(self.df["p-value"], bins=edges)
        bins = edges[1:]
        self.pval_hist = pd.DataFrame({"bin": bins, "count": counts})

        if table_filter is not None:
            self.df = self.df[self.df.apply(table_filter, axis=1)]

    def get_pval_cdf_points(self, pvals):
        idx = self.pval_cdf["p-value"].searchsorted(pvals, side="right") - 1
        d = self.pval_cdf.iloc[idx]
        d.reset_index(drop=True, inplace=True)
        return d

    def plot_pvals(self, control_targets=None, mode="hide"):
        """
        Plot the gene ranking in form of their p-values as CDF plot.
        """

        fdr_idx = self.df["fdr"].searchsorted([0.01, 0.05, 0.25], side="right")
        fdrs = pd.DataFrame(self.get_pval_cdf_points(self.df.iloc[fdr_idx]["log10-p-value"]))
        fdrs.loc[:, "label"] = self.df.iloc[fdr_idx]["fdr"].apply("{:.0%} FDR".format).values

        top5 = self.df.index
        if control_targets:
            if mode == "hide":
                valid = self.df["target"].apply(lambda target: target not in
                                                control_targets)
                top5 = top5[valid]
            elif mode == "show-only":
                valid = self.df["target"].apply(lambda target: target in
                                                control_targets)
                top5 = top5[valid]
        top5targets = top5[:5]
        top5 = pd.DataFrame(self.get_pval_cdf_points(self.df.ix[top5targets, "log10-p-value"]))
        top5.loc[:, "target"] = top5targets.values

        plt = templates.get_template("plots/pvals.json").render(
            pvals=self.pval_cdf.to_json(orient="records"),
            highlight=top5.to_json(orient="records"),
            fdrs=fdrs.to_json(orient="records"))
        return plt

    def get_pvals_highlight_targets(self, highlight_targets):
        data = pd.DataFrame(self.get_pval_cdf_points(self.df.ix[highlight_targets, "log10-p-value"]))
        data.loc[:, "target"] = highlight_targets
        return data

    def plot_pval_hist(self):
        return templates.get_template("plots/pval_hist.json").render(
            hist=self.pval_hist.to_json(orient="records"))

    def ids(self, fdr):
        valid = self.df["fdr"] <= fdr
        return set(self.df.ix[valid, "target"])

    def __len__(self):
        return self.df.shape[0]
