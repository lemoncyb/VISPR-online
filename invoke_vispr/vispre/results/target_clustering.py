from operator import itemgetter
import math

import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.cluster.vq import kmeans2

from ..results.common import templates


class TargetClustering:
    def __init__(self, target_results, topn=500):
        self.df = pd.DataFrame()
        for condition, selections in sorted(target_results.items(), key=itemgetter(0)):
            targets = pd.concat([targets[:]["score"]
                                 for targets in selections.values()])
            self.df[condition] = targets[~targets.index.duplicated()]
        self.conditions = self.df.columns
        mean = self.df.abs().mean(axis=1)
        mean.sort_values(ascending=False, inplace=True)
        self.df = self.df.ix[mean[:topn].index]
        self.linkage = linkage(self.df, method="ward", metric="euclidean")

    def plot_clustering(self, k):
        """Plot k flat clusters."""
        clustering = pd.DataFrame(fcluster(self.linkage, k,
                                           criterion="maxclust"),
                                  columns=["cluster"])
        clustering.index = self.df.index[clustering.index]
        data = pd.merge(self.df, clustering, left_index=True, right_index=True)
        data.sort_values("cluster", inplace=True)
        data["target"] = data.index

        conditions = []
        for condition in self.df.columns:
            d = data.loc[:, ["cluster", "target", condition]]
            d.columns = ["cluster", "target", "beta"]
            d.loc[:, "condition"] = condition
            conditions.append(d)
        data = pd.concat(conditions)
        betamax = math.ceil(data["beta"].abs().max())
        plt = templates.get_template("plots/target_clustering.json").render(
            data=data.to_json(orient="records"),
            beta_min=-betamax,
            beta_max=betamax)
        return plt
