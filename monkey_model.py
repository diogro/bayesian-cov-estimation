import pandas as pd
import numpy as np
import pymc as pm
import noise_control as nc
import dendropy
import operator

t = dendropy.Tree.get_from_path("./small.nwm.tree.nw", "newick")
num_leafs = len(t.leaf_nodes())
num_traits = 39
effects = ['SUB', 'SEX']

data = pd.read_csv("./monkey.data.csv")
raw_matrices = pd.read_csv("./monkey.matrices.csv")
with open('monkey.matrices.labels.txt') as f:
    monkey_labels = f.read().splitlines()

# Lendo matrizes ML pra todo mundo, junto com tamanhos amostrais


def make_symetric(matrix):
    matrix = np.tril(matrix) + np.tril(matrix, k=-1).transpose()
    return matrix

node_matrices = {monkey_labels[0]: np.array(raw_matrices.ix[0:num_traits-1, ])}
node_sample_size = {monkey_labels[0]: sum(data['species'] == monkey_labels[1])}
for i in range(1, len(monkey_labels)):
    new_matrix = np.array(raw_matrices.ix[i*num_traits:(((i+1)*num_traits)-1), :])
    node_matrices[monkey_labels[i]] = make_symetric(new_matrix)
    node_sample_size[monkey_labels[i]] = sum(data['species'] == monkey_labels[i])

# Tirando quem nao esta na filogenia e trocando os keys

for i in range(len(monkey_labels)):
    if t.find_node_with_taxon_label(monkey_labels[i]):
        new_key = str(t.find_node_with_taxon_label(monkey_labels[i]))
        node_sample_size[new_key] = node_sample_size.pop(monkey_labels[i])
        if node_sample_size[new_key] < 60:
            node_matrices[new_key] = make_symetric(nc.noise_control(node_matrices.pop(monkey_labels[i])))
        else:
            node_matrices[new_key] = node_matrices.pop(monkey_labels[i])
    else:
        node_matrices.pop(monkey_labels[i])
        node_sample_size.pop(monkey_labels[i])

# Funcao que recebe uma lista de filhos e calcula a matriz media pro parent
# node


def matrix_mean(child_labels):
    new_matrix = node_sample_size[str(child_labels[0])]*node_matrices[str(child_labels[0])]
    sample = node_sample_size[str(child_labels[0])]
    for i in range(1, len(child_labels)):
        new_matrix = new_matrix +\
            node_sample_size[str(child_labels[i])] * node_matrices[str(child_labels[i])]
        sample = sample + node_sample_size[str(child_labels[i])]
    new_matrix = make_symetric(new_matrix/sample)
    return new_matrix, sample

# Calculando as matrizes e tamanhos amostrais para todos os nodes

for n in t.postorder_node_iter():
    if str(n) not in node_matrices:
        node_matrices[str(n)], node_sample_size[str(n)] = matrix_mean(n.child_nodes())

# Agora comeca o PyMC

root = t.seed_node

theta = [pm.MvNormalCov('theta_0',
                        #mu=np.array(data.ix[:, 0:num_traits].mean()),
                        mu=np.zeros(num_traits),
                        C=np.eye(num_traits)*10.,
                        value=np.zeros(num_traits))]

sigma = [pm.WishartCov('sigma_0',
                       n=num_traits+1,
                       C=np.eye(num_traits)*10.,
                       value=node_matrices[str(root)])]

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
                               value=node_matrices[str(n)]))

    tree_idx[str(n)] = len(theta) - 1
    i = i + 1


def mk_fixed_effects(effects):
    factor_effects = {}
    for n in t.leaf_nodes():
        factor_effects[str(n)] = mk(n, effects[0], effects[1:])
    return factor_effects


def mk(s, e0, es):
    filtered = set(data.ix[data['species'] == str(s.taxon), e0])
    if len(es) < 1:
        return {k: {} for k in filtered}
    else:
        return {k: mk(s, es[0], es[1:] if len(es) > 1 else []) for k in filtered}

data_list = []
#data_sim_list = []

effects_tree = mk_fixed_effects(effects)


def tree_flattening(tree):
    tree_list = []

    def tf(t, lpart, lresult):
        if not t.items():
            lresult.append(lpart)
        else:
            map(lambda k: tf(t[k], lpart + [k], lresult), t.keys())

    tf(tree, [], tree_list)
    return tree_list

flat_effects_tree = tree_flattening(effects_tree)


def mk_node(species, node_name, node, parent_idx, effects, path, has_siblings=False):
    paths = reduce(lambda x, y: "{}__{}".format(x, y).replace(' ', '_'), path)

    if has_siblings:
        theta.append(pm.MvNormalCov('theta_{}'.format(paths),
                                    mu=theta[parent_idx],
                                    C=np.eye(num_traits),
                                    value=np.zeros(num_traits)))
        sigma.append(pm.WishartCov('sigma_{}'.format(paths),
                                   n=num_traits+1,
                                   C=sigma[parent_idx],
                                   value=node_matrices[species]))

        parent_idx = len(theta) - 1

    if not node.items():
        obs_data = np.array(data.ix[(data['species'] == str(n.taxon)) &
   	        reduce(operator.iand,
                map(lambda s: data[s[0]] == s[1],
                    zip(effects, path[1:]))), 0:num_traits])

        data_list.append(pm.MvNormalCov('data_{}'.format(paths),
                                            mu=theta[parent_idx],
                                            C=sigma[parent_idx],
                                            value=obs_data,
                                            observed=True))

        # Slow method to simulate n populations from posterior	
        #ds = []
        #for i in xrange(0,obs_data.shape[0]):
        #    ds.append(pm.MvNormalCov('data_{}_{}'.format(paths, i),
        #                                    mu=theta[parent_idx],
        #                                    C=sigma[parent_idx]
        #                                    ))

        #data_sim_list.append(ds) 

        return

    has_siblings = len(node.keys()) > 1
    for k in node.keys():
        mk_node(species, k, node[k], parent_idx, effects, path + [k], has_siblings)


for n in t.leaf_nodes():
    leaf_idx = tree_idx[str(n)]
    path = [str(n.taxon)]
    has_siblings = effects_tree[str(n)] and len(effects_tree[str(n)]) > 1
    mk_node(str(n), str(n), effects_tree[str(n)], leaf_idx, effects, path, has_siblings)
