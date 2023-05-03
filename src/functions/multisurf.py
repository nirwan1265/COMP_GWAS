# Packages
import pandas as pd
from skrebate import MultiSURF
import numpy as np
import os
os.chdir("/share/maize/ntanduk/multisurf")

# Getting the data
df = pd.read_csv("Final_filtered_data.csv")
df = df.iloc[:,:]

# Response
pheno_y = df["OlsenP"].values
y = np.log1p(pheno_y)

# Predictors
cols_to_drop = ['Database ID', 'Country', 'BIOME_NAME', 'Groups', 'Original Database', 'Converted from', 'Averaged or not', 'Depth', 'Countryf', 'Scode', 'Continent', 'ResidualsClass']
df = df.drop(columns=cols_to_drop, axis = 1)
X = df.values.astype(np.float64)

# Scaling
scaler = StandardScaler()
X = scaler.fit_transform(X)

# Model fit
ms = MultiSURF(n_jobs=-1)
ms.fit(X, y)

# Extracting the features
feature_scores = pd.DataFrame({"Feature": df.columns, "Score": ms.feature_importances_})
feature_scores.sort_values("Score", ascending=False, inplace=True)

# Saving as csv
feature_scores.to_csv("feature_scores.csv", index=False)

