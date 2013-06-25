import pymc as pm
import numpy as np

# Criando parametros conhecidos
original_sigma = np.array([[1, 0.5], [0.5, 1]])
original_theta = [1, -1]

# Simulando dados com media theta e covariancia sigma
original_theta_1 = np.random.multivariate_normal(original_theta, original_sigma)
original_theta_2 = np.random.multivariate_normal(original_theta, original_sigma)

#original_sigma_1 = pm.distributions.wishart_cov_like(original_sigma,
                                               #C=original_sigma,
                                               #n=3)
#original_sigma_2 = pm.distributions.WishartCov(original_sigma,
                                               #C=original_sigma,
                                               #n=3)

data_1 = np.random.multivariate_normal(original_theta_1, original_sigma, 100)
data_2 = np.random.multivariate_normal(original_theta_2, original_sigma, 100)

# Definindo priors, Wishart pra sigma e normal pra theta
sigma = pm.WishartCov('sigma', n=3, C=np.eye(2), value=np.eye(2))
sigma_1 = pm.WishartCov('sigma_1', n=3, C=sigma, value=np.eye(2))
sigma_2 = pm.WishartCov('sigma_2', n=3, C=sigma, value=np.eye(2))

theta = pm.MvNormalCov('theta', mu=[0.,0.], C=np.eye(2), value=[0.,0.])
theta_1 = pm.MvNormalCov('theta_1', mu=theta, C=sigma, value=[0.,0.])
theta_2 = pm.MvNormalCov('theta_2', mu=theta, C=sigma, value=[0.,0.])

# Verossimilhanca gaussiana com media theta e covariancia sigma
x_1 = pm.MvNormalCov('x_1', theta_1, sigma_1, value=data_1, observed=True)
x_2 = pm.MvNormalCov('x_2', theta_2, sigma_2, value=data_2, observed=True)
