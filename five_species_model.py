import pandas as pd
import numpy as np
import pymc as pm

dados = pd.read_csv("./dados5sp.csv")

# Filogenia: (B, ((C, E),(A, D)))

# Hiper parametros do no (B, (CEAD))

theta_B_CEAD = pm.MvNormalCov('theta_B_CEAD',
                              mu=np.array(dados.ix[:,0:4].mean()),
                              C=np.eye(4)*10.)

sigma_B_CEAD = pm.WishartCov('sigma_B_CEAD',
                             n=5,
                             C=np.eye(4)*10.)

# Ramos do no (B, (CEAD))

theta_B = pm.MvNormalCov('theta_B',
                         mu=theta_B_CEAD,
                         C=sigma_B_CEAD)

sigma_B = pm.WishartCov('sigma_B',
                        n=5,
                        C=sigma_B_CEAD)

theta_CEAD = pm.MvNormalCov('theta_CEAD',
                            mu=theta_B_CEAD,
                            C=sigma_B_CEAD)

sigma_CEAD = pm.WishartCov('sigma_CEAD',
                           n=5,
                           C=sigma_B_CEAD)

# Ramos do no CEAD

theta_CE = pm.MvNormalCov('theta_CE',
                          mu=theta_CEAD,
                          C=sigma_CEAD)

sigma_CE = pm.WishartCov('sigma_CE',
                         n=5,
                         C=sigma_CEAD)

theta_AD = pm.MvNormalCov('theta_AD',
                          mu=theta_CEAD,
                          C=sigma_CEAD)

sigma_AD = pm.WishartCov('sigma_AD',
                         n=5,
                         C=sigma_CEAD)

# Ramos do no CE

theta_C = pm.MvNormalCov('theta_C',
                         mu=theta_CE,
                         C=sigma_CE)

sigma_C = pm.WishartCov('sigma_C',
                        n=5,
                        C=sigma_CE)

theta_E = pm.MvNormalCov('theta_E',
                         mu=theta_CE,
                         C=sigma_CE)

sigma_E = pm.WishartCov('sigma_E',
                        n=5,
                        C=sigma_CE)

# Ramos do no AD

theta_A = pm.MvNormalCov('theta_A',
                         mu=theta_AD,
                         C=sigma_AD,
                         value=np.zeros(4))

sigma_A = pm.WishartCov('sigma_A',
                        n=5,
                        C=sigma_AD)

theta_D = pm.MvNormalCov('theta_D',
                         mu=theta_AD,
                         C=sigma_AD)

sigma_D = pm.WishartCov('sigma_D',
                        n=5,
                        C=sigma_AD)

# Verossimilhancas dos terminais
data_A = pm.MvNormalCov('data_A',
                        theta_A,
                        sigma_A,
                        value=np.array(dados.ix[dados['especies']=='A', 0:4]),
                        observed=True)
data_B = pm.MvNormalCov('data_B',
                        theta_B,
                        sigma_B,
                        value=np.array(dados.ix[dados['especies']=='B', 0:4]),
                        observed=True)
data_C = pm.MvNormalCov('data_C',
                        theta_C,
                        sigma_C,
                        value=np.array(dados.ix[dados['especies']=='C', 0:4]),
                        observed=True)
data_D = pm.MvNormalCov('data_D',
                        theta_D,
                        sigma_D,
                        value=np.array(dados.ix[dados['especies']=='D', 0:4]),
                        observed=True)
data_E = pm.MvNormalCov('data_E',
                        theta_E,
                        sigma_E,
                        value=np.array(dados.ix[dados['especies']=='E', 0:4]),
                        observed=True)
