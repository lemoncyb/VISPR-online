# coding: utf-8
from __future__ import absolute_import, division, print_function

__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster, Liu lab"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import os

import pandas as pd

from ..results import target
from ..results import rna
from ..results import fastqc
from ..results import mapstats
from ..results import target_overlap
from ..results import target_clustering
from ..common import VisprError


class Screens(object):
    def __init__(self):
        self.screens = {}

    def add(self, config, parentdir="."):
        screen = config["experiment"]
        self.screens[screen] = Screen(config, parentdir=parentdir)

    def __iter__(self):
        return map(self.__getitem__, sorted(self.screens.keys())) # [cuiyb] for Python 3.5
        #return iter(map(self.__getitem__, sorted(self.screens.keys()))) # [cuiyb] for Python 2.7

    def __getitem__(self, screen):
        return self.screens[screen]

    def _overlap_targets(self, fdr=0.05, items=None):
        assert items is not None

        def label(screen, condition, selection):
            if condition != "default":
                return screen, condition, selection
            else:
                return screen, selection

        return {
            "\n".join(label(screen, condition, selection)):
            self.screens[screen].targets[condition][selection].ids(fdr)
            for screen, condition, selection in items
        }

    def plot_overlap_venn(self, fdr=0.05, items=None):
        plt = target_overlap.plot_overlap_venn(
            **self._overlap_targets(fdr=fdr,
                                    items=items))
        return plt

    def overlap(self, fdr=0.05, items=None):
        return target_overlap.overlap(
            *self._overlap_targets(fdr=fdr,
                                   items=items).values())


class Screen(object):
    def __init__(self, config, parentdir="."):
        def get_path(relpath):
            if relpath is None:
                return None
            if relpath.startswith("http"):
                return relpath
            return os.path.join(parentdir, relpath)

        self.name = config["experiment"]

        if self.name=='mle':
            self.targets, self.is_mle = parse_target_results(
                get_path(config["targets"]["results"]))
        elif self.name=='bagel':
            self.targets, self.is_mle = parse_target_results_bagel(
                get_path(config["targets"]["results"]))
        elif self.name=='jacks':
            self.targets, self.is_mle = parse_target_results_jacks(
                get_path(config["targets"]["results"]),get_path(config["sgrnas"]["counts"]),get_path(config["sgrnas"]["results"]))                
        else:
            self.targets, self.is_mle = parse_target_results(
                get_path(config["targets"]["results"]))

        self.is_genes = config["targets"].get("genes", False)
        self.species = config["species"].upper()
        self.save = config["save"]
        if "session" in config:
            self.session = config["session"]
        #self.assembly = config["assembly"]  # used for IGV

        self.rnas = rna.Results(
            get_path(config["sgrnas"]["counts"]),
            info=get_path(config["sgrnas"].get("annotation", None)),
            posterior_results=get_path(config["sgrnas"].get("results", None)))
        self.mapstats = None
        if "mapstats" in config["sgrnas"]:
            self.mapstats = mapstats.Results(
                get_path(config["sgrnas"]["mapstats"]))

        self.fastqc = None
        if "fastqc" in config:
            self.fastqc = fastqc.Results(**{
                sample: map(get_path, paths)
                for sample, paths in config["fastqc"].items()
            })

        self.control_targets = set()
        if "controls" in config["targets"]:
            try:
                self.control_targets = set(
                    pd.read_table(get_path(config["targets"]["controls"]),
                                  header=None,
                                  squeeze=True,
                                  na_filter=False,
                                  low_memory=False))
            except (Exception, BaseException) as e:
                raise VisprError(
                    "Failed to parse control targets "
                    "(targets->controls in config): {}".format(e))

        self.target_clustering = None
        if self.is_mle:
            # provide clustering on score
            self.target_clustering = target_clustering.TargetClustering(self.targets)

    def targets(self, positive=True):
        return self.pos_targets if positive else self.neg_targets


def parse_target_results(path,
                         selections=["negative selection",
                                     "positive selection"]):
    results = pd.read_table(path, na_filter=False, low_memory=False)
    paths = [col.split("|") for col in results.columns]
    is_mle = "beta" in [path[1] for path in paths if len(path) > 1]

    if is_mle:
        # MLE format
        def get_results(condition, selection):
            res = results[[results.columns[0]] +
                          "{cond}|beta {cond}|p-value {cond}|fdr".format(
                              cond=condition).split()]
            res.columns = ["target", "score", "p-value", "fdr"]
            if selection == "negative selection":
                table_filter = lambda row: row["score"] <= 0
            else:
                table_filter = lambda row: row["score"] >= 0
            return target.Results(res.copy(), table_filter=table_filter)

        conditions = [path[0] for path in paths if len(path) > 1]
 
        targets = {
            condition: {
                selection: get_results(condition, selection)
                for selection in selections
            }
            for condition in conditions
        }

        return targets, True
    else:
        # RRA format
        def get_results(selection):
            res = results[[results.columns[0]] +
                          "{sel}|score {sel}|p-value {sel}|fdr".format(
                              sel=selection[:3]).split()]
            res.columns = ["target", "score", "p-value", "fdr"]
            return target.Results(res.copy())

        targets = {
            "default":
            {selection: get_results(selection)
             for selection in selections}
        }
        return targets, False

def parse_target_results_bagel(path,
                         selections=["foldchange"]):
    results = pd.read_table(path, na_filter=False, low_memory=False)

    def get_results(selection):
        sum = [results.columns[0]]
        for i in range(0, results.shape[1]):
                sum = sum + [results.columns[i]]
        res = results[sum]

        tmpCols = ["target"]
        for col in results.columns.values.tolist():
            tmpCols.append(col)
        # if results.shape[1]>2:
        #     for s in ["bf%d"%i for i in range(1, results.shape[1])]:
        #         ls.append(s)
        res.columns = tmpCols
        return target.Results(res.copy())

    targets = {
        "default":
        {selection: get_results(selection)
         for selection in selections}
    }
    return targets, False

def parse_target_results_jacks(path,path_counts,path_grna,
                         selections=["genescore","foldchange","grna"]):
    results = pd.read_table(path, na_filter=False, low_memory=False)

    # add logfoldchange
    results_foldchange = pd.read_table(path_counts, na_filter=False, low_memory=False)
   
    # add grna
    results_grna = pd.read_table(path_grna, na_filter=False, low_memory=False)

    def get_results(selection):
        if selection == "genescore":
            sum = [results.columns[0]]
            for i in range(0, results.shape[1]):
                sum = sum + [results.columns[i]]
            res = results[sum]

            tmpCols = ["target"]
            for col in results.columns.values.tolist():
                tmpCols.append(col)
            res.columns = tmpCols    
        elif selection == "foldchange":
            sum = [results_foldchange.columns[0]]
            for i in range(0, results_foldchange.shape[1]):
                sum = sum + [results_foldchange.columns[i]]
            res = results_foldchange[sum]

            tmpCols = ["target"]
            for col in results_foldchange.columns.values.tolist():
                tmpCols.append(col)
            res.columns = tmpCols  
        else:
            sum = [results_grna.columns[0]]
            for i in range(0, results_grna.shape[1]):
                sum = sum + [results_grna.columns[i]]
            res = results_grna[sum]
            
            tmpCols = ["target"]
            for col in results_grna.columns.values.tolist():
                tmpCols.append(col)
            res.columns = tmpCols
        return target.Results(res.copy())

    targets = {
        "default":
        {selection: get_results(selection)
         for selection in selections}
    }
    return targets, False