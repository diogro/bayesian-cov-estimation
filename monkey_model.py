import pandas as pd
import numpy as np
import pymc as pm
import dendropy

data = pd.read_csv("./monkey.data.csv")

t = dendropy.Tree.get_from_path("./monkey.tree.nw", "newick")
num_leafs = len(t.leaf_nodes())
num_traits = 39

root = t.seed_node

theta = [pm.MvNormalCov('theta_0',
                        mu=np.array(data.ix[:,0:num_traits].mean()),
                        C=np.eye(num_traits)*10.,
                        value=np.zeros(num_traits))]

sigma = [pm.WishartCov('sigma_0',
                       n=num_traits+1,
                       C=np.eye(num_traits)*10.,
                       value=np.eye(num_traits))]

tree_idx = {str(root): 0}

i = 1
for n in t.nodes()[1:]:
    parent_idx = tree_idx[str(n.parent_node)]

    theta.append(pm.MvNormalCov('theta_{}'.format(str(i)),
                                mu=theta[parent_idx],
                                C=sigma[parent_idx],
                                value=np.zeros(num_traits)))

    sigma.append(pm.WishartCov('sigma_{}'.format(str(i)),
                               n=num_traits+1,
                               C=sigma[parent_idx],
                               value=np.eye(num_traits)))

    tree_idx[str(n)] = len(theta) - 1
    i = i + 1

data_list = []
for n in t.leaf_nodes():
    leaf_idx = tree_idx[str(n)]
    data_list.append(pm.MvNormalCov('data_{}'.format(n.taxon),
                                    mu=theta[leaf_idx],
                                    C=sigma[leaf_idx],
                                    value=np.array(data.ix[data['species'] == str(n.taxon), 0:num_traits]),
                                    observed=True))
