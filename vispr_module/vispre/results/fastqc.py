# coding: utf-8
from __future__ import absolute_import, division, print_function

__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster, Liu lab"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


import json
import re
from io import StringIO

import pandas as pd
from flask import render_template


from ..results.common import templates
from ..common import VisprError


class Results(object):
    def __init__(self, **fastqc_data):
        self.gc_content = []
        self.base_quality = []
        self.seq_quality = []
        self.n_samples = len(fastqc_data)
        for sample, paths in fastqc_data.items():
            for path in paths:
                filename = parse_fastqc_filename(path)
                for name, data in parse_fastqc_data(path):
                    if name == "Per sequence GC content":
                        data.columns = ["gc", "density"]
                        data["density"] /= data["density"].sum()
                        self.gc_content.append(data)
                    elif name == "Per base sequence quality":
                        data.columns = ["base", "mean", "median", "lo", "hi", "q10", "q90"]
                        self.base_quality.append(data)
                    elif name == "Per sequence quality scores":
                        data.columns = ["qual", "density"]
                        data["density"] /= data["density"].sum()
                        self.seq_quality.append(data)
                    data["sample"] = sample
                    data["filename"] = filename
        if self.gc_content:
            self.gc_content = pd.concat(self.gc_content)
        else:
            self.gc_content = None
        if self.base_quality:
            self.base_quality = pd.concat(self.base_quality)
        else:
            self.base_quality = None
        if self.seq_quality:
            self.seq_quality = pd.concat(self.seq_quality)
        else:
            self.seq_quality = None

    def plot_gc_content(self):
        plt = ""
        if self.gc_content is not None:
            plt = templates.get_template("plots/gc_content.json").render(
                data=self.gc_content.to_json(orient="records"),
                show_legend=self.n_samples <= 20)
        return plt

    def plot_base_quality(self):
        plt = ""
        if self.base_quality is not None:
            plt = templates.get_template(
                "plots/base_quality.json").render(
                data=self.base_quality.to_json(orient="records"),
                show_legend=self.n_samples <= 20)
        return plt

    def plot_seq_quality(self):
        plt = ""
        if self.seq_quality is not None:
            plt = templates.get_template(
                "plots/seq_quality.json").render(
                data=self.seq_quality.to_json(orient="records"),
                show_legend=self.n_samples <= 20)
        return plt


def parse_fastqc_filename(path, pattern=re.compile(r"Filename\t(?P<name>.*)")):
    with open(path) as f:
        match = pattern.search(f.read())
        if match:
            return match.group("name")
    raise IOError("Invalid FastQC format.")


def parse_fastqc_data(path, pattern=re.compile(r"\n\>\>(?P<name>[\w ]+)\t(pass|fail)\n\#(?P<data>.+?)\n\>\>END_MODULE", flags=re.MULTILINE | re.DOTALL)):
    with open(path) as f:
        f = f.read()
        for match in pattern.finditer(f):
            d = match.group("data")
            try:
                d = StringIO(d)
            except TypeError:
                d = StringIO(unicode(d))
            try:
                data = pd.read_table(d, low_memory=False)
            except (Exception, BaseException) as e:
                raise VisprError(
                    "Failed to parse FastQC results: {}".format(e))
            yield match.group("name"), data
