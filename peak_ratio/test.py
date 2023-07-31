# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 19:08:43 2023

@author: worms
"""

# Import necessary libraries
import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols
import matplotlib.pyplot as plt
import seaborn as sns

# Load the data
data = pd.read_csv("/mnt/data/baranyi_res_2.csv")

# Run a two-way ANOVA for the 'lag' response variable
model_lag = ols('lag ~ C(enrichment) + C(induction) + C(enrichment):C(induction)', data=data).fit()
anova_lag = sm.stats.anova_lm(model_lag, typ=2)

# Run a two-way ANOVA for the 'mumax' response variable
model_mumax = ols('mumax ~ C(enrichment) + C(induction) + C(enrichment):C(induction)', data=data).fit()
anova_mumax = sm.stats.anova_lm(model_mumax, typ=2)

# Print ANOVA results
print(anova_lag)
print(anova_mumax)

# Change the labels for better interpretation
data['induction'] = data['induction'].replace({'glu': 'Glucose', 'ara': 'Arabinose'})

# Create boxplots
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

sns.boxplot(x='enrichment', y='lag', data=data, ax=axes[0, 0])
axes[0, 0].set_title('Effect of Enrichment on Lag')

sns.boxplot(x='enrichment', y='mumax', data=data, ax=axes[0, 1])
axes[0, 1].set_title('Effect of Enrichment on Mumax')

sns.boxplot(x='induction', y='lag', data=data, ax=axes[1, 0])
axes[1, 0].set_title('Effect of Induction on Lag')

sns.boxplot(x='induction', y='mumax', data=data, ax=axes[1, 1])
axes[1, 1].set_title('Effect of Induction on Mumax')

plt.tight_layout()
plt.show()
