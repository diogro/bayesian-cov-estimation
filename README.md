Estimando matrizes de covariância com modelos hierárquicos
==========================================================

Estratégia geral
----------------

Utilizar métodos hierárquicos para estimar matrizes de covariância e
médias fenotípicas de forma filogeneticamente estruturada.

Modelo
------

Dada um filogenia e um conjunto de medidas quantitativas, podemos escrever um modelo hierárquico para estimar médias e covariâncias entre as medidas utilizando aproximações gaussianas.

Passos do modelo:

1. Notação:
    + $x_ij$ matriz de medidas, i populações e j caracteres.
    + $\theta_i$ médias das populações
    + $\sigma_i$ matriz de covariância das populações
    + $\Theta_k$ hiper parâmetros, médias nos nós da filogenia
    + $\Sigma_k$ hiper parâmetros, matriz de covariâncias nos nós da filogenia
2. Verossimilhança das populações
    + gaussiana básica: $p(z|\theta, \sigma) = \sum_{ij} exp\left(-\frac{1}{2}*(x_{ij} - \theta_i)\sigma_i^1(x_{ij}-\theta_i)\right)$
3. Priors hierárquicos
