# Pegasos Kernel SVM Implementation

> 🧠 **Machine Learning** · 📊 **Binary Classification** · 📐 **Support Vector Machines** · 🔬 **Kernel Methods**

## 📌 Project Overview

This project explores the **Pegasos algorithm** for training a binary Support Vector Machine (SVM) using two kernel functions:

- **Linear kernel** — for learning linear decision boundaries.
- **Radial Basis Function (RBF) kernel** — for learning nonlinear decision boundaries.

The original Pegasos algorithm uses stochastic or mini-batch updates. This implementation instead evaluates the entire training dataset at each iteration, producing a deterministic, full-batch variant when the data and initialization remain fixed.

The project focuses on understanding the mathematical foundations of kernel SVMs, implementing the optimization procedure, and examining how different kernels affect classification performance.

### 🎥 Video Tutorials

- [Video Tutorial 1](https://www.youtube.com/watch?v=Q7vT0--5VII)
- [Video Tutorial 2](https://www.youtube.com/watch?v=N_RQj4OL1mg)
- [Video Tutorial 3](https://www.youtube.com/watch?v=OKFMZQyDROI)

## 🎯 Main Goal

Implement and understand a kernel-based SVM training procedure inspired by Pegasos, including margin violations, hinge loss, regularization, and kernel-based predictions.

## 🎯 Project Objectives

- Understand the mathematical principles behind binary SVM classification.
- Implement linear and RBF kernel functions.
- Construct a kernel matrix for efficient training-set calculations.
- Understand the role of the hinge loss and margin condition.
- Explore regularization and the Pegasos learning-rate schedule.
- Implement deterministic full-batch updates.
- Compare linear and nonlinear decision boundaries.
- Evaluate classification performance using appropriate metrics.

## 🧭 Methodology

| Step | Description |
|---|---|
| 1. Prepare the data | Format the input features and encode labels as −1 and +1. |
| 2. Select a kernel | Choose a linear or RBF kernel. |
| 3. Construct the kernel matrix | Calculate pairwise similarities between training samples. |
| 4. Initialize parameters | Initialize the dual coefficients and bias. |
| 5. Calculate decision scores | Evaluate the current model on all training samples. |
| 6. Identify margin violations | Find samples for which the signed decision score is below 1. |
| 7. Update parameters | Apply the specified regularization and hinge-loss updates. |
| 8. Evaluate the model | Measure classification performance on held-out data. |

## 🧮 1. Binary Classification and Decision Function

For binary classification, the target labels are encoded as:

\[
y_i \in \{-1,+1\}
\]

The kernel-based decision function is:

\[
f(x)=\sum_{j=1}^{N}\alpha_j y_j K(x_j,x)+b
\]

where:

- \(x_j\) is a training sample.
- \(y_j\) is its class label.
- \(\alpha_j\) is its nonnegative dual coefficient.
- \(K(x_j,x)\) measures similarity between two samples.
- \(b\) is the bias term.
- \(N\) is the number of training samples.

The predicted class is:

\[
\hat y=\operatorname{sign}(f(x))
\]

A positive decision score predicts class \(+1\), while a negative score predicts class \(-1\). A score of exactly zero lies on the decision boundary, so the implementation should define how such ties are handled.

## 📐 2. SVM Margin Condition

The signed decision score is:

\[
y_i f(x_i)
\]

For a correctly classified sample outside or on the margin:

\[
y_i f(x_i)\geq 1
\]

A sample violates the margin when:

\[
y_i f(x_i)<1
\]

Margin violations include both misclassified samples and correctly classified samples that lie inside the margin.

## 📉 3. Hinge Loss

The hinge loss for a single sample is:

\[
L_i=\max(0,1-y_i f(x_i))
\]

Its interpretation is straightforward:

- If \(y_i f(x_i)\geq 1\), the hinge loss is zero.
- If \(y_i f(x_i)<1\), the hinge loss is positive.

The hinge loss penalizes samples that are misclassified or insufficiently separated from the decision boundary.

A standard linear SVM objective combines regularization and average hinge loss:

\[
\min_w\frac{\lambda}{2}\|w\|^2+
\frac{1}{N}\sum_{i=1}^{N}
\max(0,1-y_i(w^\top x_i+b))
\]

Here, \(\lambda>0\) controls the strength of regularization. A kernelized implementation must express the corresponding regularization and update rules consistently in terms of its coefficients.

## ⏱️ 4. Pegasos Learning Rate

The Pegasos learning-rate schedule is:

\[
\eta_t=\frac{1}{\lambda t}
\]

where:

- \(\eta_t\) is the learning rate at iteration \(t\).
- \(\lambda\) is the regularization parameter.
- \(t\) is the iteration number, starting at 1.

The learning rate decreases as training progresses, allowing the optimization process to make progressively smaller updates.

**Important distinction:** using this learning-rate schedule does not, by itself, make an algorithm a standard Pegasos implementation. The parameterization, regularization, sample-selection strategy, and update equations must also be consistent with the intended objective.

## 🔑 5. The Kernel Trick

The kernel trick allows an SVM to learn nonlinear decision boundaries without explicitly constructing a high-dimensional feature representation.

Instead of directly calculating a weight vector in feature space, the model uses kernel evaluations:

\[
f(x)=\sum_{j=1}^{N}\alpha_j y_j K(x_j,x)+b
\]

The coefficients determine the contribution of the training samples to the decision function. In a standard dual SVM formulation, the coefficients are constrained to be nonnegative, with additional constraints depending on the formulation.

## 📏 6. Linear Kernel

The linear kernel calculates the dot product between two samples:

\[
K(x_i,x_j)=x_i^\top x_j
\]

For a feature matrix `X` with shape `(n_samples, n_features)`, the kernel matrix can be calculated using NumPy:

```python
K = X @ X.T
```

The resulting matrix has shape `(n_samples, n_samples)`.

**Characteristics:**

- Produces a linear decision boundary in the original feature space.
- Is relatively inexpensive to calculate.
- Is a useful baseline for comparing kernel methods.
- Does not require an RBF-specific parameter such as \(\gamma\).

## 🌐 7. RBF Kernel

The Radial Basis Function kernel is defined as:

\[
K(x_i,x_j)=
\exp\left(-\gamma\|x_i-x_j\|^2\right)
\]

where \(\gamma>0\) controls how quickly similarity decreases with squared distance.

For the usual positive values of \(\gamma\):

- Identical samples have a kernel value of 1.
- Nearby samples have kernel values close to 1.
- More distant samples have kernel values closer to 0.

**Characteristics:**

- Can represent nonlinear decision boundaries.
- Can capture relationships that a linear kernel cannot.
- Requires the choice of a suitable \(\gamma\).
- Can be sensitive to feature scaling and hyperparameter selection.

Feature standardization should generally be considered before applying the RBF kernel.

## 🧩 8. Kernel Matrix

For \(N\) training samples, the kernel matrix is:

\[
K_{ij}=K(x_i,x_j)
\]

or, in matrix form:

\[
K=
\begin{bmatrix}
K(x_1,x_1)&K(x_1,x_2)&\cdots\\
K(x_2,x_1)&K(x_2,x_2)&\cdots\\
\vdots&\vdots&\ddots
\end{bmatrix}
\]

The matrix has shape \(N\times N\).

For a symmetric kernel, such as the linear and RBF kernels, \(K_{ij}=K_{ji}\).

Precomputing the matrix avoids repeatedly evaluating the same training-sample similarities during optimization. However, storing the full matrix requires \(O(N^2)\) memory, which can become expensive for large datasets.

## ⚙️ 9. Deterministic Full-Batch Training

In this project, every iteration evaluates all training samples rather than selecting a random sample or mini-batch.

For training sample \(x_i\), the decision score is:

\[
f_i=\sum_{j=1}^{N}\alpha_jy_jK(x_j,x_i)+b
\]

The signed margin is:

\[
m_i=y_i f_i
\]

A margin-violation mask can be computed as:

```python
violations = y * scores < 1
```

The algorithm then uses the violating samples to calculate the hinge-loss contribution to the update.

A mathematically consistent implementation must specify:

1. How regularization is applied to the model parameters.
2. How the violating samples contribute to the coefficient updates.
3. How the class labels enter those updates.
4. Whether the hinge-loss contribution is averaged over all training samples.
5. How the bias is updated or optimized.
6. Whether any coefficient constraints or projection steps are required.

These details are essential because a full-batch update is not automatically equivalent to the original stochastic Pegasos update.

### Understanding the Important Code

#### Kernel contribution

```python
scores = K @ (alpha * y) + bias
```

This calculates the decision scores for all training samples, assuming `K[i, j] = K(x_i, x_j)`.

In mathematical notation:

\[
f(x_i)=\sum_{j=1}^{N}K(x_i,x_j)\alpha_jy_j+b
\]

For a symmetric kernel, this is equivalent to using \(K(x_j,x_i)\).

#### Margin check

```python
violations = y * scores < 1
```

This produces a Boolean array identifying samples that violate the SVM margin.

For example:

```text
y      = [ 1, -1,  1]
scores = [ 2, 0.2, 0.5]
```

The signed margins are:

```text
y * scores = [2, -0.2, 0.5]
```

The resulting violation mask is:

```text
[False, True, True]
```

The second and third samples violate the margin.

#### Regularization

A shrinkage step may appear in an update as:

```python
alpha *= (1 - eta * lambda_param)
```

However, **this line should not be assumed to be the correct regularization step for every kernelized Pegasos formulation**. The appropriate regularization update depends on how the model is parameterized and how the primal or dual objective is derived.

#### Updating violating samples

An update such as:

```python
alpha[violations] += eta / n_samples
```

is not a complete, generally valid kernel SVM update on its own. In particular, the update must account for the labels and the chosen objective, and the bias requires a defined update rule.

Use the update equations derived for the specific algorithm rather than treating these lines as a complete implementation.

## 🔄 10. Algorithm Flow

```text
Training data X, y
        |
        v
Validate features and encode labels
        |
        v
Choose linear or RBF kernel
        |
        v
Construct the training kernel matrix K
        |
        v
Initialize coefficients and bias
        |
        v
For each iteration
        |
        +--> Calculate learning rate
        |
        +--> Compute decision scores
        |
        +--> Calculate signed margins
        |
        +--> Identify margin violations
        |
        +--> Apply the derived regularization update
        |
        +--> Apply the derived hinge-loss update
        |
        +--> Update or optimize the bias
        |
        v
Return the trained model
        |
        v
Evaluate on held-out data
```

## ⚖️ 11. Linear vs. RBF Kernel

| Feature | Linear Kernel | RBF Kernel |
|---|---|---|
| Kernel function | \(K(x,z)=x^\top z\) | \(K(x,z)=e^{-\gamma\|x-z\|^2}\) |
| Decision boundary | Linear | Potentially nonlinear |
| Main advantage | Simplicity and efficiency | Flexible nonlinear modeling |
| Additional parameter | No \(\gamma\) | Requires \(\gamma\) |
| Feature scaling | Often beneficial | Particularly important |
| Computational cost | Usually lower | Can be higher |
| Suitable use case | Approximately linear patterns | More complex nonlinear patterns |

The better kernel depends on the dataset, feature representation, regularization, and hyperparameter choices. It should be selected using validation data rather than assumed in advance.

## 📊 12. Model Evaluation

The implementation can be evaluated using a held-out test set and suitable binary-classification metrics, including:

- **Accuracy:** overall proportion of correct predictions.
- **Precision:** proportion of positive predictions that are correct.
- **Recall:** proportion of actual positive samples correctly identified.
- **F1-score:** harmonic mean of precision and recall.
- **Confusion matrix:** breakdown of correct and incorrect predictions.

If the dataset is imbalanced, precision, recall, F1-score, and the confusion matrix can provide more useful insight than accuracy alone.

Hyperparameters such as \(\lambda\), \(\gamma\), and the number of iterations should be selected using training and validation data. The test set should be reserved for final evaluation.

## 🛠️ Technologies & Tools

| Technology | Purpose |
|---|---|
| Python | Implementation |
| NumPy | Vectorized calculations and kernel matrices |
| Matplotlib | Visualization, if used |
| scikit-learn | Data preparation and evaluation utilities, if used |
| Jupyter Notebook | Interactive experimentation, if used |

Keep this list aligned with the libraries actually used in the project.

## 📁 Project Structure

```text
Pegasos-Kernel-SVM/
│
├── README.md
├── project.ipynb
└── requirements.txt
```

This is a suggested structure; adjust it to match the actual repository files.

## 🚀 Expected Outcomes

By completing the project, the implementation should demonstrate:

- Understanding of the SVM margin and hinge loss.
- An implementation of linear and RBF kernel calculations.
- Efficient computation of decision scores using a kernel matrix.
- A clear distinction between stochastic Pegasos and full-batch optimization.
- An evaluation of classification performance on held-out data.
- An understanding of the trade-offs between linear and nonlinear kernels.

## 🔬 Scientific Motivation

Kernel SVMs provide a useful framework for studying how similarity functions influence classification. The linear kernel operates directly on the original feature representation, while the RBF kernel enables nonlinear decision boundaries through pairwise similarities.

Implementing these methods from their mathematical foundations helps connect optimization, regularization, margin-based learning, and kernel methods. It also highlights the importance of deriving parameter updates from a well-defined objective rather than relying on individual update statements in isolation.

## 📌 Project Focus

**Data preparation → Kernel construction → Margin-based optimization → Model prediction → Performance evaluation → Linear/RBF comparison**

The central prediction rule is:

\[
\boxed{
f(x)=\sum_{i=1}^{N}\alpha_i y_i K(x_i,x)+b
}
\]

Classification is then based on the sign of the decision score, with an explicit tie-handling rule for scores equal to zero.
