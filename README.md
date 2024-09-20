# PySmilarityMol

## Introduction

Welcome to the **PySmilarityMol** repository, where we explore molecular similarity and clustering. This project demonstrates methods for comparing compounds using **Tanimoto similarity**, **Euclidean distance**, and molecular descriptors (e.g., molecular weight, logP). The project also showcases **hierarchical clustering** to group similar compounds, visualized via dendrograms. This is particularly useful in fields like drug discovery and chemical analysis.

## Step 1: Calculating and Visualizing Compound Similarity

This script calculates and visualizes the similarity between compounds from a CSV file containing SMILES strings.

### Key Steps:
1. **SMILES Input**: Load the compounds' SMILES strings.
2. **Morgan Fingerprints**: Generate fingerprints using RDKit.
3. **Tanimoto Similarity**: Calculate the Tanimoto similarity matrix.
4. **Euclidean Distance**: Use molecular weight as a descriptor to calculate the Euclidean distance.
5. **Heatmap Visualization**: Visualize both matrices as heatmaps.

```python
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import pairwise_distances

# Load SMILES data
df = pd.read_csv("smiles.csv")

# Generate Morgan Fingerprints
fps = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smiles), 2) for smiles in df['SMILES']]

# Calculate Tanimoto Similarity
tanimoto_matrix = [[DataStructs.TanimotoSimilarity(fps[i], fps[j]) for j in range(len(fps))] for i in range(len(fps))]
df_tanimoto = pd.DataFrame(tanimoto_matrix)

# Calculate Euclidean Distance
df['MolWt'] = df['SMILES'].apply(lambda x: Descriptors.MolWt(Chem.MolFromSmiles(x)))
euclidean_matrix = pairwise_distances(df[['MolWt']])
df_euclidean = pd.DataFrame(euclidean_matrix)

# Plot Heatmaps
sns.heatmap(df_tanimoto, cmap="YlGnBu", annot=True)
plt.title("Tanimoto Similarity")
plt.show()

sns.heatmap(df_euclidean, cmap="YlGnBu", annot=True)
plt.title("Euclidean Distance")
plt.show()

Step 2: Hierarchical Clustering and Dendrogram

Using the Tanimoto similarity matrix, we apply hierarchical clustering to identify groups of similar compounds. The clustering result is visualized as a dendrogram.

python
Always show details

from scipy.cluster import hierarchy

# Perform hierarchical clustering
Z = hierarchy.linkage(df_tanimoto, method='ward')

# Visualize Dendrogram
plt.figure(figsize=(25, 10))
dn = hierarchy.dendrogram(Z)
plt.title("Hierarchical Clustering (Tanimoto Similarity)")
plt.show()

Step 3: Euclidean Distance Clustering

Similarly, we perform hierarchical clustering using the Euclidean distance matrix to cluster compounds based on molecular weight.

python
Always show details

# Perform clustering based on Euclidean distance
Z = hierarchy.linkage(df_euclidean, method='ward')

# Dendrogram visualization
plt.figure(figsize=(25, 10))
dn = hierarchy.dendrogram(Z)
plt.title("Hierarchical Clustering (Euclidean Distance)")
plt.show()

License

This repository is licensed under the MIT License. Please see the LICENSE file for more details.

MIT License

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
Contact

    Email: rikgangulybioinfo@gmail.com
    Computational Biology Laboratory, North-Eastern Hill University, Shillong, India """

