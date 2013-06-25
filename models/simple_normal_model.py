import pymc as pm
import numpy as np

# Criando parametros conhecidos
original_sigma = np.array([[1, 0.5], [0.5, 1]])
original_theta = [1, -1]

# Simulando dados com media theta e covariancia sigma
data = np.random.multivariate_normal(original_theta, original_sigma, 100)

# Definindo priors, Wishart pra sigma e normal pra theta
sigma = pm.WishartCov('sigma', n=3, C=np.eye(2), value=np.eye(2))

theta = pm.MvNormalCov('theta', mu=[0.,0.], C=np.eye(2), value=[0.,0.])

# Verossimilhanca gaussiana com media theta e covariancia sigma
x = pm.MvNormalCov('x', theta, sigma, value=data, observed=True)
