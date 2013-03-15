import pandas as pd
import numpy as np
import pymc as pm
import dendropy as dendro

data = pd.read_csv("./dados5sp.csv")

t = dendropy.Tree.get_from_string("(a,(b,(c,d)))", "newick")
num_leafs = len(t.leaf_nodes())
num_traits = 4

root = t.seed_node

theta = [pm.MVNormalCov('theta_0',
			mu=np.array(data.ix[:,0:num_traits].mean()),
			C=np.eye(num_traits)*10.)]

sigma = [pm.WishartCov('sigma_0',
			n=num_leafs,
			C=np.eye(num_traits)*10.)]

tree_idx = {str(root): 0}

for n in t.nodes()[1:]:
    parent_idx = tree_idx(str(n.parent_node))

    theta.append('theta_{}'.format(parent_idx),
		 mu=theta[parent_idx],
		 C=sigma[parent_idx])

    sigma.append('sigma_{}'.format(parent_idx),
		 n=num_leafs,
		 C=sigma[parent_idx])
 
    tree_idx[str(n)] = len(theta) - 1

data = []
i = 0  
for n in t.leaf_nodes():
    leaf_idx = tree_idx[str(n)]
    data.append(pm.MvNormalCov('data_{}'.format(i),
			      theta[leaf_idx],
			      sigma[leaf_idx],
		              value=np.array(data.ix[data['species']==n.taxon.label, 0:num_traits),
			      observed=True))
    i += 1
    
