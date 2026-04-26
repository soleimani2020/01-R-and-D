# AdaBoost From Scratch in Python 🚀

## Overview

This repository contains a **pure Python/NumPy implementation of AdaBoost (Adaptive Boosting)** – no scikit-learn or other high‑level ML libraries.  

AdaBoost is a foundational ensemble method that combines multiple **weak learners** (here, decision stumps) into a **strong classifier** by iteratively focusing on previously misclassified samples.

### Core Concepts Covered

- Weak learner optimization (decision stumps)  
- Weighted sample training  
- Weighted error minimization  
- Learner importance (α) calculation  
- Iterative boosting rounds  
- Ensemble classifier construction  

---

## Key Features ✨

- ✅ Pure Python + NumPy (no ML libraries)  
- 🌲 Decision stump weak learners  
- 🎯 Threshold & polarity optimization  
- ⚖️ Weighted error minimization  
- 📈 Adaptive sample weight updates  
- 🔁 Configurable number of boosting iterations  
- 🧠 Clean, educational code for learning  

---

## How AdaBoost Works ⚙️

For each boosting round \( t = 1 \dots T \):

### 1️⃣ Initialize Sample Weights (first round only)

All \( N \) samples start with equal importance:

\[
w_i^{(1)} = \frac{1}{N}, \quad i = 1, \dots, N
\]

### 2️⃣ Train a Weak Learner (Decision Stump)

Find the best stump (feature \( j \), threshold \( \theta \), polarity \( p \in \{-1, 1\} \)) that minimizes the **weighted error**.

### 3️⃣ Compute Weighted Error

\[
\text{error}_t = \sum_{i=1}^{N} w_i^{(t)} \cdot \mathbf{1}\big(y_i \neq h_t(\mathbf{x}_i)\big)
\]

Where:
- \( w_i^{(t)} \) = weight of sample \( i \) at round \( t \)
- \( y_i \) = true label (\( \pm 1 \))
- \( h_t(\mathbf{x}_i) \) = prediction of weak learner \( t \)
- \( \mathbf{1}(\cdot) \) = indicator function

### 4️⃣ Compute Learner Importance (Alpha)

\[
\alpha_t = \frac{1}{2} \ln\left(\frac{1 - \text{error}_t}{\text{error}_t}\right)
\]

- Lower error → higher \( \alpha \) → more influence in the final vote.

### 5️⃣ Update Sample Weights

\[
w_i^{(t+1)} = w_i^{(t)} \, \exp\!\big(-\alpha_t \, y_i \, h_t(\mathbf{x}_i)\big)
\]

- Correctly classified samples get lower weight; misclassified get higher weight.

### 6️⃣ Normalize Weights

\[
w_i^{(t+1)} \leftarrow \frac{w_i^{(t+1)}}{\sum_{j=1}^{N} w_j^{(t+1)}}
\]

### 7️⃣ Repeat

Continue for \( T \) boosting rounds (`n_clf`).

---

## Final Strong Classifier 🏆

The ensemble prediction for a new sample \( \mathbf{x} \) is:

\[
H(\mathbf{x}) = \text{sign}\!\left( \sum_{t=1}^{T} \alpha_t \, h_t(\mathbf{x}) \right)
\]

Where `sign` returns `+1` or `-1`.

---

## Example Dataset 📊

```python
import numpy as np

X = np.array([
    [1, 2],
    [2, 3],
    [3, 4],
    [4, 5]
])
y = np.array([1, 1, -1, -1])
