import pandas as pd
import numpy as np
import pymc as pm
import dendropy

import os, sys, inspect
# realpath() with make your script run, even if you symlink it :)
cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../utils")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)
import noise_control as nc

t = dendropy.Tree.get_from_path("../trees/nwm.genus.tree.nw", "newick")
num_leafs = len(t.leaf_nodes())
num_traits = 39

data = pd.read_csv("../dados/mean.center.residuals.NWM.csv")

matrices = data.groupby('genus').apply(lambda x: x.cov())
genus = pd.unique(data['genus'])
# Lendo matrizes ML pra todo mundo, junto com tamanhos amostrais

def make_symetric(matrix):
    matrix = np.tril(matrix) + np.tril(matrix, k=-1).transpose()
    return matrix

node_matrices = {}
node_sample_size = {}
for g in genus:
    new_matrix = np.array(matrices.T[g])
    #print new_matrix.shape
    #node_matrices[g] = make_symetric(new_matrix)
    node_matrices[g] = new_matrix
    node_sample_size[g] = sum(data['genus'] == g)

# Tirando quem nao esta na filogenia e trocando os keys

node_means = {}

node_name = lambda n: str(n.label or n.taxon or n)

for g in genus:
    if t.find_node_with_taxon_label(g):
        new_key = node_name(t.find_node_with_taxon_label(g))
        node_means[new_key] = np.array(data.ix[data['genus'] == str(g), 0:num_traits]).mean(0)
        node_sample_size[new_key] = node_sample_size.pop(g)
        if node_sample_size[new_key] < num_traits + 2:
            node_matrices[new_key] = make_symetric(nc.noise_control(node_matrices.pop(g)))
        else:
            node_matrices[new_key] = node_matrices.pop(g)
    else:
        node_matrices.pop(g)
        node_sample_size.pop(g)

# Funcao que recebe uma lista de filhos e calcula a matriz media pro parent
# node

def matrix_mean(child_labels):
    new_matrix = np.zeros((num_traits, num_traits))
    sample = 0
    new_mean = np.zeros(num_traits)
    for child in child_labels:
        new_matrix = new_matrix +\
            node_sample_size[node_name(child)] * node_matrices[node_name(child)]
        sample = sample + node_sample_size[node_name(child)]
        new_mean = new_mean + node_means[node_name(child)]
    new_matrix = new_matrix/sample
    new_mean = new_mean/len(child_labels)
    return new_matrix, sample, new_mean

# Calculando as matrizes e tamanhos amostrais para todos os nodes

for n in t.postorder_node_iter():
    if node_name(n) not in node_matrices:
        node_matrices[node_name(n)], node_sample_size[node_name(n)], node_means[node_name(n)] = matrix_mean(n.child_nodes())

# Agora comeca o PyMC

root = t.seed_node

theta = {node_name(root): pm.MvNormalCov('theta_0',
                                          mu=np.array(data.ix[:, 0:num_traits].mean()),
                                          #value=node_means[str(root)],
                                          value=np.zeros(num_traits),
                                          #mu=np.zeros(num_traits),
                                          C=np.eye(num_traits)*100.)}

sigma = {node_name(root): pm.WishartCov('sigma_0',
                                         value=node_matrices[node_name(root)],
                                         #value=np.eye(num_traits),
                                         n=num_traits+1,
                                         C = node_matrices[node_name(root)])}
#C=np.eye(num_traits)*100.)}

#var_factors = {}
betas = {}

for n in t.nodes()[1:]:
    parent_idx = node_name(n.parent_node)

    #var_factors[str(i)] = pm.Uniform('var_factor_{}'.format(str(i)), lower=0, upper=1000)

    betas[node_name(n)] = pm.MvNormalCov('betas_{}'.format(node_name(n)),
                                         value=np.zeros(num_traits),
                                         mu=np.zeros(num_traits),
                                         C=np.eye(num_traits)*100.)

    theta[node_name(n)] = pm.MvNormalCov('theta_{}'.format(node_name(n)),
                                         value=node_means[node_name(n)],
                                         #value=np.zeros(num_traits),
                                         #mu=theta[parent_idx],
                                         mu=theta[parent_idx] + betas[node_name(n)],
                                         C=sigma[parent_idx])
    #C=sigma[parent_idx]*var_factors[node_name(n)])

    sigma[node_name(n)] = pm.WishartCov('sigma_{}'.format(node_name(n)),
                                        value=node_matrices[node_name(n)],
                                        #value=np.eye(num_traits),
                                        n=num_traits+1,
                                        C=sigma[parent_idx])

data_list = []
for n in t.leaf_nodes():
    leaf_idx = node_name(n)
    data_list.append(pm.MvNormalCov('data_{}'.format(n.taxon),
                                    mu=theta[leaf_idx],
                                    C=sigma[leaf_idx],
                                    value=np.array(data.ix[data['genus'] == str(n.taxon), 0:num_traits]),
                                    observed=True))
