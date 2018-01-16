# coding: utf-8
from __future__ import absolute_import, division, print_function

__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster, Liu lab"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import re, json, os
import string
import random, uuid
import shutil

import numpy as np
from flask import Flask, render_template, request, session, abort, flash, redirect, url_for
import yaml
from werkzeug.utils import secure_filename
from pymongo import MongoClient

from vispr_screen import __version__

from vispr_screen import app
from vispr_screen.results import Screens
from vispr_screen.common import VisprError

app.jinja_env.trim_blocks = True
app.jinja_env.lstrip_blocks = True
app.secret_key = ''.join(
        random.choice(string.ascii_uppercase + string.digits)
        for _ in range(30))

client = MongoClient('localhost', 27017)    #Configure the connection to the database
db = client.vispr    #Select the database
sequence = db.sequence #Select the collection
annotation = db.annotation

d3cate10 = ["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf"]
google20c = ["#3366CC","#DC3912","#FF9900","#109618","#990099","#3B3EAC","#0099C6","#DD4477","#66AA00","#B82E2E",
             "#316395","#994499","#22AA99","#AAAA11","#6633CC","#E67300","#8B0707","#329262","#5574A6","#3B3EAC"]

with open(os.path.join(os.path.dirname(__file__), "captions.yaml")) as f:
    CAPTIONS = yaml.load(f)

UPLOAD_FOLDER = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))),'tmp')


def init_server(*configs):
    app.screens = Screens()
    print("Loading data.")
    for path in configs:
        with open(path) as f:
            config = yaml.load(f)
        try:
            app.screens.add(config, parentdir=os.path.dirname(path))
        except KeyError as e:
            raise VisprError(
                "Syntax error in config file {}. Missing key {}.".format(path,
                                                                         e))
    app.secret_key = ''.join(
        random.choice(string.ascii_uppercase + string.digits)
        for _ in range(30))


@app.route("/ping")
def ping():
    return ""


@app.route("/", methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        species = request.form.get('species')
        file_gene = request.files["gene"]
        file_count = request.files["count"]
        if "sgrna" in request.files:
            file_sgrna = request.files["sgrna"]
        else:
            file_sgrna = ""
        if "library" in request.files:
            file_lib = request.files["library"]
        else:
            file_lib = ""
        save_session = request.form.get("save")

        tmp_dir = str(uuid.uuid4())
        try:
            os.makedirs(os.path.join(UPLOAD_FOLDER, tmp_dir))
        except KeyError as e:
            "Error when creating tmp directory.".format(e)

        for file in [file_gene, file_count, file_sgrna, file_lib]:
            if file:
                filename = file.filename
                file.save(os.path.join(UPLOAD_FOLDER, tmp_dir, filename))

        vispr_config = {
            "experiment": 'mle',
            "species": species,
            "targets": {
                "results": file_gene.filename,
                "genes": "true"
            },
            "sgrnas": {
                "counts": file_count.filename
            }
        }
        if file_sgrna:
            vispr_config["sgrnas"]["results"] = file_sgrna.filename
        if file_lib:
            vispr_config["sgrnas"]["annotation"] = file_lib.filename
        if save_session:
            vispr_config["save"] = True
            vispr_config["session"] = tmp_dir
        else:
            vispr_config["save"] = False
        config_file = os.path.join(UPLOAD_FOLDER, tmp_dir, 'vispr.yaml')
        with open(config_file, "w") as f:
            yaml.dump(vispr_config, f, default_flow_style=False)

        init_server(config_file)
        if not save_session:
            shutil.rmtree(os.path.join(UPLOAD_FOLDER, tmp_dir))  # delete uploaded files
        screen = next(iter(app.screens))
        #condition = next(iter(screen.targets))
        #return  redirect(url_for('targets', screen=screen.name, condition=condition, selection='positive selection'))
        return  redirect(url_for('target_clustering', screen=screen.name))

    elif hasattr(app, 'screens'):
        screen = next(iter(app.screens))
        return render_template("index.html",
                               screens=app.screens,
                               screen=screen,
                               version=__version__)
    else:
        return render_template("index.html",
                           screens=None,
                           screen=None,
                           version=__version__)


@app.route("/<screen>")
def index_screen(screen):
    try:
        screen = app.screens[screen]
    except KeyError:
        abort(404)
    return render_template("index.html",
                           screens=app.screens,
                           screen=screen,
                           version=__version__)


@app.route("/targets/<screen>/<condition>/<selection>")
def targets(screen, condition, selection):
    screen = app.screens[screen]
    table_args = request.query_string.decode()

    gorilla = screen.species in GORILLA_SPECIES and screen.is_genes
    targets = ""
    background = ""
    if gorilla:
        overlap_args = get_overlap_args()
        targets = screen.targets[condition][selection][:]["target"]
        if overlap_args:
            overlap = app.screens.overlap(*overlap_args)
            filter = targets.apply(lambda target: target in overlap)
            background = targets[~filter]
            targets = targets[filter]
            background = background.to_csv(None, index=False)
        targets = targets.to_csv(None, index=False)

    return render_template("targets.html",
                           captions=CAPTIONS,
                           screens=app.screens,
                           selection=selection,
                           condition=condition,
                           screen=screen,
                           control_targets=screen.control_targets,
                           hide_control_targets=session.get(
                               "control_targets_mode", "hide") == "hide",
                           table_args=table_args,
                           samples=screen.rnas.samples,
                           has_rna_info=screen.rnas.info is not None,
                           gorilla=gorilla,
                           gorilla_targets=targets,
                           gorilla_background=background,
                           gorilla_mode="hg" if background else "mhg",
                           version=__version__)


@app.route("/target-clustering/<screen>")
def target_clustering(screen):
    screen = app.screens[screen]
    return render_template("target_clustering.html",
                           screens=app.screens,
                           screen=screen,
                           captions=CAPTIONS,
                           version=__version__)


@app.route("/qc/<screen>")
def qc(screen):
    screen = app.screens[screen]
    return render_template("qc.html",
                           captions=CAPTIONS,
                           screens=app.screens,
                           screen=screen,
                           fastqc=screen.fastqc is not None,
                           mapstats=screen.mapstats is not None,
                           version=__version__)


@app.route("/compare/<screen>")
def compare(screen):
    screen = app.screens[screen]
    return render_template("compare.html",
                           captions=CAPTIONS,
                           screens=app.screens,
                           screen=screen,
                           version=__version__)


def tbl_targets(screen, condition, selection,
                offset=None,
                perpage=None,
                add_locus=False):
    screen = app.screens[screen]

    # sort and slice records
    records = screen.targets[condition][selection][:]
    total_count = records.shape[0]
    filter_count = total_count

    filter = np.ones(total_count, dtype=np.bool)

    # restrict to overlap
    overlap_args = get_overlap_args()
    if overlap_args:
        overlap = app.screens.overlap(*overlap_args)
        filter &= records["target"].apply(lambda target: target in overlap)

    # searching
    search = get_search()
    if search:
        filter &= records["target"].str.contains(search)

    control_targets_mode = session.get("control_targets_mode", "hide")
    control_targets = screen.control_targets
    if control_targets_mode == "hide":
        filter &= records["target"].apply(
            lambda target: target not in control_targets)
    elif control_targets_mode == "show-only":
        filter &= records["target"].apply(
            lambda target: target in control_targets)

    # filtering
    if not np.all(filter):
        if np.any(filter):
            records = records[filter]
            filter_count = records.shape[0]
        else:
            return render_template("dyntable.json",
                                   records="[]",
                                   filter_count=0,
                                   total_count=total_count)

    # sorting
    columns, ascending = get_sorting()
    if columns:
        records = records.sort(columns, ascending=ascending)
    if offset is not None:
        records = records[offset:offset + perpage]

    # formatting
    def fmt_col(col):
        if col.dtype == np.float64:
            return col.apply("{:.2g}".format)
        return col

    records = records.apply(fmt_col)

    if add_locus and screen.rnas.info is not None:

        def get_locus(target):
            locus = screen.rnas.target_locus(target)
            if locus is None:
                return ""
            return "{}:{}-{}".format(*locus)

        records["locus"] = list(map(get_locus, records["target"]))

    return records, filter_count, total_count


@app.route("/tbl/targets/json/<screen>/<condition>/<selection>",
           methods=["GET"])
def tbl_targets_json(screen, condition, selection):
    offset = int(request.args.get("offset", 0))
    perpage = int(request.args.get("perPage", 20))
    records, filter_count, total_count = tbl_targets(screen, condition,
                                                     selection,
                                                     offset=offset,
                                                     perpage=perpage,
                                                     add_locus=False)
    return render_template("dyntable.json",
                           records=records.to_json(orient="records",
                                                   double_precision=15),
                           filter_count=filter_count,
                           total_count=total_count)


@app.route("/tbl/targets/txt/<screen>/<condition>/<selection>",
           methods=["GET"])
def tbl_targets_txt(screen, condition, selection):
    records, _, _ = tbl_targets(screen, condition, selection)
    return records[["target", "score", "p-value", "fdr"]].to_csv(sep="\t",
                                                                 index=False)


@app.route("/tbl/pvals_highlight/<screen>/<condition>/<selection>/<targets>")
def tbl_pvals_highlight(screen, condition, selection, targets):
    screen = app.screens[screen]
    targets = targets.split("|")
    records = screen.targets[condition][selection].get_pvals_highlight_targets(
        targets)
    return records.to_json(orient="records")


def gene_seq(target):
    return sequence.find_one({"name":target})


def gene_anno(target):
    return annotation.find({"gene_name":target})


def add_locus(screen, target):
    seq = gene_seq(target)
    if seq is None:
        return False
    trans = gene_anno(target)
    if trans is None:
        return False

    def cal_locus(tran):
        return {"id": tran["transcript_id"],
                "exons": list(map(lambda x,y,z: {"x":x, "y":y, "description":z, "color":"#22AA99"},
                                  tran["exon_start"], tran["exon_stop"], tran["exon_id"]))
                }

    rnas_locus = screen.rnas.rnas_locus(target)

    return {"gene_id": seq["id"],
            "gene_start": seq["start"],
            "gene_stop": seq["stop"],
            "gene_seq": seq["seq"],
            "trans": list(map(cal_locus, trans)),
            "rnas": list(map(lambda x,y,z,c: {"x":str(x), "y":str(y), "description":z, "color":c},
                             rnas_locus["start"], rnas_locus["stop"], rnas_locus["rna"], google20c))
            }


@app.route("/tbl/locus/<screen>/<target>")
def tbl_locus(screen, target):
    screen = app.screens[screen]
    locus = add_locus(screen, target)
    return json.dumps(locus)


@app.route("/tbl/rnas/<screen>/<target>")
def tbl_rnas(screen, target):
    screen = app.screens[screen]
    table = screen.rnas.by_target(target)
    table.sort_index(inplace=True)
    table["color"] = google20c[0:len(table.index)]
    return table.to_json(orient="records")


@app.route("/igv/session/<screen>.xml")
def igv_session(screen):
    screen = app.screens[screen]
    return render_template("igv/session.xml", screen=screen)


@app.route("/igv/track/rnas/<screen>.gct")
def igv_track_rnas(screen):
    screen = app.screens[screen]
    return screen.rnas.track()


@app.route("/plt/pvals/<screen>/<condition>/<selection>")
def plt_pvals(screen, condition, selection):
    screen = app.screens[screen]
    plt = screen.targets[condition][selection].plot_pvals(
        screen.control_targets,
        mode=session.get("control_targets_mode", "hide"))
    return plt


@app.route("/plt/pvalhist/<screen>/<condition>/<selection>")
def plt_pval_hist(screen, condition, selection):
    screen = app.screens[screen]
    plt = screen.targets[condition][selection].plot_pval_hist()
    return plt


@app.route("/plt/normalization/<screen>")
def plt_normalization(screen):
    screen = app.screens[screen]
    plt = screen.rnas.plot_normalization()
    return plt


@app.route("/plt/pca/<screen>/<int:x>/<int:y>/<int:legend>")
def plt_pca(screen, x, y, legend):
    screen = app.screens[screen]
    plt = screen.rnas.plot_pca(comp_x=x, comp_y=y, legend=legend == 1, )
    return plt


@app.route("/plt/correlation/<screen>")
def plt_correlation(screen):
    screen = app.screens[screen]
    plt = screen.rnas.plot_correlation()
    return plt


@app.route("/plt/gc_content/<screen>")
def plt_gc_content(screen):
    screen = app.screens[screen]
    plt = screen.fastqc.plot_gc_content()
    return plt


@app.route("/plt/base_quality/<screen>")
def plt_base_quality(screen):
    screen = app.screens[screen]
    plt = screen.fastqc.plot_base_quality()
    return plt


@app.route("/plt/seq_quality/<screen>")
def plt_seq_quality(screen):
    screen = app.screens[screen]
    plt = screen.fastqc.plot_seq_quality()
    return plt


@app.route("/plt/mapstats/<screen>")
def plt_mapstats(screen):
    screen = app.screens[screen]
    plt = screen.mapstats.plot_mapstats()
    return plt


@app.route("/plt/zerocounts/<screen>")
def plt_zerocounts(screen):
    screen = app.screens[screen]
    plt = screen.mapstats.plot_zerocounts()
    return plt


@app.route("/plt/gini_index/<screen>")
def plt_gini_index(screen):
    screen = app.screens[screen]
    plt = screen.mapstats.plot_gini_index()
    return plt


@app.route("/plt/readcount_cdf/<screen>")
def plt_readcounts(screen):
    screen = app.screens[screen]
    plt = screen.rnas.plot_readcount_cdf()
    return plt


@app.route("/plt/target-clustering/<screen>/<int:k>")
def plt_target_clustering(screen, k):
    screen = app.screens[screen]
    plt = screen.target_clustering.plot_clustering(k)
    return plt


@app.route("/plt/overlap_chord")
def plt_overlap_chord():
    return app.screens.plot_overlap_chord(*get_overlap_args())


@app.route("/plt/overlap_venn")
def plt_overlap_venn():
    return app.screens.plot_overlap_venn(*get_overlap_args())


@app.route("/set/control_targets_mode/<mode>")
def set_control_targets_mode(mode):
    session["control_targets_mode"] = mode
    return ""


def get_overlap_args():
    if "fdr" not in request.values and "overlap-items" not in request.form:
        return None

    fdr = float(request.values.get("fdr", 0.25))
    items = list(map(lambda item: item.split("|"),
                     request.values.getlist("overlap-items")))
    return fdr, items


def get_sorting(pattern=re.compile("sorts\[(?P<col>.+)\]")):
    cols, ascending = [], []
    for arg, val in request.args.items():
        m = pattern.match(arg)
        if m:
            cols.append(m.group("col"))
            ascending.append(int(val) == 1)
    return cols, ascending


def get_search(pattern=re.compile("search\[(?P<target>.+)\]")):
    return request.args.get("queries[search]", None)


def get_targets(screen, selection):
    return screen.targets(selection == "positive")


GORILLA_SPECIES = (
    "HOMO_SAPIENS ARABIDOPSIS_THALIANA SACCHAROMYCES_CEREVISIAE "
    "CAENORHABDITIS_ELEGANS DROSOPHILA_MELANOGASTER DANIO_RERIO "
    "MUS_MUSCULUS RATTUS_NORVEGICUS".split())
