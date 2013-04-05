import pandas as pd
import numpy as np
import pymc as pm
import dendropy

data = pd.read_csv("./monkey.data.csv")
raw_matrices = pd.read_csv("./monkey.matrices.csv")
with open('monkey.matrices.labels.txt') as f:
    monkey_labels = f.read().splitlines()

num_traits = 39
node_matrices = {monkey_labels[0]: np.array(raw_matrices.ix[0:num_traits-1,])}
node_sample_size = {monkey_labels[0]: sum(data['species'] == monkey_labels[1])}
for i in range(1, len(monkey_labels)):
    node_matrices[monkey_labels[i]] = np.array(raw_matrices.ix[i*num_traits:(((i+1)*num_traits)-1),:])
    node_sample_size[monkey_labels[i]] = sum(data['species'] == monkey_labels[i])

t = dendropy.Tree.get_from_path("./small.nwm.tree.nw", "newick")
num_leafs = len(t.leaf_nodes())

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

sub_effects = {}
for n in t.leaf_nodes():
    sub_list = list(data.ix[data['species']==str(n.taxon),'SUB'])
    sub_effects[str(n)] = [sub_list[1]]
    for sub in sub_list[1:]:
        if not (any(sub in s for s in sub_effects[str(n)])):
            sub_effects[str(n)].append(sub)

sex_effects = {}
for n in t.leaf_nodes():
    sex_list = list(data.ix[data['species']==str(n.taxon),'SEX'])
    sex_effects[str(n)] = [sex_list[1]]
    for sex in sex_list[1:]:
        if not (any(sex in s for s in sex_effects[str(n)])):
            sex_effects[str(n)].append(sex)

data_list = []
for n in t.leaf_nodes():
    leaf_idx = tree_idx[str(n)]
    for sub in sub_effects[str(n)]:
        for sex in sex_effects[str(n)]:
            theta.append(pm.MvNormalCov('theta_{}_{}_{}'.format(n.taxon, str(sub), sex),
                                        mu=theta[leaf_idx],
                                        C=np.eye(num_traits),
                                        value=np.zeros(num_traits)))
            data_list.append(pm.MvNormalCov('data_{}_{}_{}'.format(n.taxon, str(sub), sex),
                                            mu=theta[len(theta)-1],
                                            C=sigma[leaf_idx],
                                            value=np.array(data.ix[(data['species'] == str(n.taxon)) & (data['SUB'] == sub) & (data['SEX'] == sex), 0:num_traits]),
                                            observed=True))
