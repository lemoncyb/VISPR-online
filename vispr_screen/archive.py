import sys
import os
import tarfile
from io import BytesIO

import yaml

from vispr_screen.common import VisprError


def archive(config, out):
    parentdir = os.path.dirname(config)

    def get_path(relpath):
        if relpath.startswith("http"):
            return relpath
        return os.path.join(parentdir, relpath)

    with open(config) as config:
        config = yaml.load(config)

    def get_tar_path(name):
        return os.path.join(config["experiment"], name)

    mode = None
    if out.endswith(".tar"):
        mode = "w"
    elif out.endswith(".tar.gz"):
        mode = "w:gz"
    elif out.endswith(".tar.bz2"):
        mode = "w:bz2"

    with tarfile.open(out, mode) as tar:
        try:
            if "fastqc" in config:
                new = {}
                for sample, fastqs in config["fastqc"].items():
                    new[sample] = []
                    for i, f in enumerate(fastqs):
                        out = "{}.{}.fastqc_data.txt".format(sample, i)
                        tar.add(get_path(f), get_tar_path(out))
                        new[sample].append(out)
                config["fastqc"] = new

            out = "all.count.normalized.txt"
            tar.add(get_path(config["sgrnas"]["counts"]), get_tar_path(out))
            config["sgrnas"]["counts"] = out

            if "mapstats" in config["sgrnas"]:
                out = "all.countsummary.txt"
                tar.add(get_path(config["sgrnas"]["mapstats"]), get_tar_path(out))
                config["sgrnas"]["mapstats"] = out

            if "annotation" in config["sgrnas"]:
                out = "sgnra_annotation.bed"
                tar.add(get_path(config["sgrnas"]["annotation"]), get_tar_path(out))
                config["sgrnas"]["annotation"] = out

            if "results" in config["sgrnas"]:
                out = "all.sgrna_summary.txt"
                tar.add(get_path(config["sgrnas"]["results"]), get_tar_path(out))
                config["sgrnas"]["results"] = out

            out = "all.gene_summary.txt"
            tar.add(get_path(config["targets"]["results"]), get_tar_path(out))
            config["targets"]["results"] = out

            if "controls" in config["targets"]:
                out = "all.controls.txt"
                tar.add(get_path(config["targets"]["controls"]), get_tar_path(out))
                config["targets"]["controls"] = out
        except KeyError as e:
            raise VisprError("Missing key in config file ({}). "
                             "This is most likely not a valid vispr config file. "
                             "If generated with MAGeCK-VISPR, the config file "
                             "should end with .vispr.yaml and be located under "
                             "results/ in the working directory "
                             "(e.g., results/myexperiment.vispr.yaml).".format(e))

        newconfig = yaml.dump(config, default_flow_style=False)
        tarinfo = tarfile.TarInfo(get_tar_path("vispr.yaml"))
        tarinfo.size = len(newconfig)
        tar.addfile(tarinfo, fileobj=BytesIO(newconfig.encode()))
