# coding: utf-8
from __future__ import absolute_import, division, print_function

__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster, Liu lab"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import argparse
import logging
import sys
import string
import random
import os
import shutil
import tarfile
try:
    from urllib.request import urlopen
    from urllib.error import HTTPError
except ImportError:
    from urllib2 import urlopen, HTTPError

import yaml
from appdirs import AppDirs

from vispr_module.vispr import Screens, Screen
from vispr_module.vispr.common import VisprError
#from vispr.server import app
from vispr_module import app
from vispr_module.vispr.version import __version__
from vispr_module.vispr.archive import archive as _archive

appdirs = AppDirs("VISPR", "liulab")


def init_server(configs, host="127.0.0.1", port=5000):
    app.screens = Screens()
    print("Loading data.")
    print(configs)
    #for path in configs:
    with open(configs) as f:
        config = yaml.load(f)
    try:
        app.screens.add(config, parentdir=os.path.dirname(configs))
    except KeyError as e:
        raise VisprError(
            "Syntax error in config file {}. Missing key {}.".format(configs,
                                                                         e))
    app.secret_key = ''.join(
        random.choice(string.ascii_uppercase + string.digits)
        for _ in range(30))
    logging.info("Starting server.")
    logging.info("")
    logging.info(
        "Open:  go to " + host + ":{} in your browser.".format(port))
    logging.info("Note: Safari and Internet Explorer are currently unsupported.")
    logging.info("Close: hit Ctrl-C in this terminal.")
    #app.run(host=host, port=port)


def test_server(host="127.0.0.1",port=5000, update=False):
    testdir = os.path.join(appdirs.user_cache_dir, "testdata")
    datasets = "ESC-MLE Leukemia-MLE Melanoma-MLE".split()
    for dataset in datasets:
        exists = os.path.exists(os.path.join(testdir, dataset))
        if update or not exists:
            logging.info("Downloading {} test data.".format(dataset))
            try:
                testdata = tarfile.open(fileobj=urlopen(
                    "https://bitbucket.org/liulab/"
                    "vispr/downloads/{}.testdata.tar.bz2".format(dataset)),
                                        mode="r|bz2")
            except HTTPError as e:
                logging.error("Unable to download test data. This is most likely a temporary problem with bitbucket.org.")
                logging.error(str(e))
            if exists:
                shutil.rmtree(os.path.join(testdir, dataset))
            testdata.extractall(testdir)
    init_server(*[os.path.join(testdir, dataset, "vispr.yaml")
                  for dataset in datasets],
                host=host,
                port=port)


def print_example_config():
    print(open(os.path.join(os.path.dirname(__file__),
                            "example.config.yaml")).read())


def plots(configpath, prefix):
    directory = os.path.dirname(prefix)
    if not os.path.exists(directory):
        os.makedirs(directory)

    with open(configpath) as f:
        screen = Screen(yaml.load(f), parentdir=os.path.dirname(configpath))

    def write(json, name):
        with open(prefix + name + ".vega.json", "w") as out:
            print(json, file=out)

    if screen.fastqc is not None:
        write(screen.fastqc.plot_gc_content(), "gc-content")
        write(screen.fastqc.plot_base_quality(), "base-quality")
        write(screen.fastqc.plot_seq_quality(), "seq-quality")
    if screen.mapstats is not None:
        write(screen.mapstats.plot_mapstats(), "mapped-unmapped")
        write(screen.mapstats.plot_zerocounts(), "zerocounts")
        write(screen.mapstats.plot_gini_index(), "gini-index")
    write(screen.rnas.plot_normalization(), "readcounts")
    write(screen.rnas.plot_readcount_cdf(), "readcount-cdf")
    write(screen.rnas.plot_correlation(), "correlation")
    write(screen.rnas.plot_pca(1, 2), "pca-1-2")
    write(screen.rnas.plot_pca(1, 3), "pca-1-3")
    write(screen.rnas.plot_pca(2, 3), "pca-2-3")
    for condition, results in screen.targets.items():
        for selection, results in results.items():
            pre = ".".join(([] if condition == "default" else [condition]) +
                           [selection.replace(" ", "-")])
            write(results.plot_pvals(), pre + ".p-values")
            write(results.plot_pval_hist(), pre + ".p-value-hist")


def main():
    # create arg parser
    parser = argparse.ArgumentParser(
        description=
        "An HTML5-based interactive visualization of CRISPR/Cas9 screen data.")
    parser.add_argument("--version",
                        action="store_true",
                        help="Print version info.")
    parser.add_argument("--debug",
                        action="store_true",
                        help="Print debug info.")
    subparsers = parser.add_subparsers(dest="subcommand")

    config = subparsers.add_parser(
        "config",
        description=
        "Print an example VISPR config file. Pipe the output into a file "
        "and edit it to setup a new experiment to be displayed in VISPR.")

    server = subparsers.add_parser("server",
                                   description="Start the VISPR server.")
    server.add_argument(
        "config",
        nargs="+",
        help="YAML config files. Each file points to the results of one "
        "MAGeCK run.")

    server.add_argument("--host",
                        default="127.0.0.1",
                        type=str,
                        help="Host ip location to listen for client connection.")

    server.add_argument("--port",
                        default=5000,
                        type=int,
                        help="Port to listen for client connection.")

    plot = subparsers.add_parser(
        "plot",
        description="Plot visualizations in VEGA JSON format.")
    plot.add_argument("config",
                      help="YAML config file pointing to MAGeCK results.")
    plot.add_argument("prefix",
                      help="Prefix of all resulting plots. "
                      "This can be a path to a subdirectory.")

    test = subparsers.add_parser(
        "test",
        description="Start the VISPR server with test data.")
    test.add_argument("--host",
                        default="127.0.0.1",
                        type=str,
                        help="Host ip location to listen for client connection.")
    test.add_argument("--port",
                      default=5000,
                      type=int,
                      help="Port to listen for client connection.")
    test.add_argument("--update",
                      action="store_true",
                      help="Update test data. "
                      "This will redownload all test data.")

    archive = subparsers.add_parser(
        "archive",
        description=
        "Create a compressed archive for easy distribution of a given "
        "config file with all referenced files.")
    archive.add_argument("config",
                         help="YAML config file pointing to MAGeCK results.")
    archive.add_argument("tarfile",
                         help="The tar archive to write. Has to "
                         "end with .tar, .tar.gz, or .tar.bz2")

    args = parser.parse_args()

    logging.basicConfig(format="%(message)s",
                        level=logging.DEBUG if args.debug else logging.INFO,
                        stream=sys.stderr)
    logging.getLogger('werkzeug').setLevel(logging.DEBUG if args.debug else logging.ERROR)

    try:
        if args.version:
            print(__version__)
            exit(0)
        if args.subcommand == "server":
            init_server(*args.config, host=args.host, port=args.port)
        elif args.subcommand == "test":
            test_server(port=args.port, host=args.host, update=args.update)
        elif args.subcommand == "config":
            print_example_config()
        elif args.subcommand == "plot":
            plots(args.config, args.prefix)
        elif args.subcommand == "archive":
            _archive(args.config, args.tarfile)
        else:
            parser.print_help()
            exit(1)
    except VisprError as e:
        logging.error(e)
        exit(1)
    except ImportError as e:
        logging.error("{}. Please ensure that all dependencies from "
                      "requirements.txt are installed.".format(e))
        exit(1)
    exit(0)
