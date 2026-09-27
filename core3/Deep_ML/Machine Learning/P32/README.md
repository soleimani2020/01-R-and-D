# Pegasos Kernel SVM Implementation

## Overview

Video Tutorial : https://www.youtube.com/watch?v=Q7vT0--5VII

Video Tutorial2 : https://www.youtube.com/watch?v=N_RQj4OL1mg

Video Tutorial3 : https://www.youtube.com/watch?v=OKFMZQyDROI

This task implements a **deterministic Pegasos algorithm** for training a binary Support Vector Machine (SVM) with either:

- a **linear kernel**
- an **RBF kernel**

The original Pegasos algorithm is stochastic and usually selects one random training sample at each iteration.

Here, however, we use **all training samples in every iteration**, making the algorithm deterministic.

---

## Binary Classification

The labels should be:

$$
y_i \in \{-1, +1\}
$$

For each sample $x_i$, the classifier computes:

$$
f(x_i)
=
\sum_{j=1}^{N}
\alpha_j y_j K(x_j,x_i)
+
b
$$

The predicted class is determined by:

$$
\hat y_i
=
\operatorname{sign}(f(x_i))
$$

---

# 1. SVM Margin Condition

For a correctly classified sample outside the margin:

$$
y_i f(x_i) \ge 1
$$

A sample violates the margin when:

$$
y_i f(x_i) < 1
$$

These samples contribute to the hinge-loss gradient.

---

# 2. Hinge Loss

The SVM hinge loss is:

$$
L_i
=
\max(0,1-y_i f(x_i))
$$

Therefore:

- if $y_i f(x_i) \ge 1$, the sample has no hinge loss
- if $y_i f(x_i) < 1$, the sample contributes to the update

---

# 3. Pegasos Learning Rate

At iteration $t$, Pegasos uses:

$$
\eta_t
=
\frac{1}{\lambda t}
$$

where:

- $\eta_t$ = learning rate
- $\lambda$ = regularization parameter
- $t$ = iteration number

As training continues, the learning rate becomes smaller.

---

# 4. Kernel Trick

Instead of explicitly computing a weight vector $w$, a kernel SVM represents the decision function using training samples:

$$
f(x)
=
\sum_{j=1}^{N}
\alpha_j y_j K(x_j,x)
+
b
$$

The coefficients $\alpha_j$ determine how strongly each training sample contributes to the decision boundary.

---

# 5. Linear Kernel

The linear kernel is simply the dot product:

$$
K(x_i,x_j)
=
x_i^T x_j
$$

In NumPy:

```python
K = X @ X.T
```

---

# 6. RBF Kernel

The Radial Basis Function kernel is:

$$
K(x_i,x_j)
=
\exp
\left(
-\gamma
\|x_i-x_j\|^2
\right)
$$

where $\gamma$ controls how quickly similarity decreases with distance.

A small distance gives:

$$
K(x_i,x_j) \approx 1
$$

A large distance gives:

$$
K(x_i,x_j) \approx 0
$$

---

# 7. Kernel Matrix

For $N$ training samples, we construct an $N \times N$ kernel matrix:

$$
K =
\begin{bmatrix}
K(x_1,x_1) & K(x_1,x_2) & \cdots \\
K(x_2,x_1) & K(x_2,x_2) & \cdots \\
\vdots & \vdots & \ddots
\end{bmatrix}
$$

This allows us to calculate the decision values for all samples efficiently.

---

# 8. Deterministic Pegasos Update

Unlike the original stochastic Pegasos algorithm, we evaluate **every sample at every iteration**.

First calculate:

$$
f_i
=
\sum_j
\alpha_j y_j K(x_j,x_i)
+
b
$$

Then check the margin:

$$
y_i f_i < 1
$$

If the condition is true, sample $i$ violates the margin.

The algorithm updates the coefficients based on all margin-violating samples.

---



# Understanding the Important Lines

## Kernel contribution

```python
scores = K @ (alpha * y) + bias
```

This implements:

$$
f(x_i)
=
\sum_j
K(x_i,x_j)
\alpha_j y_j
+
b
$$

---

## Margin Check

```python
violations = y * scores < 1
```

This checks:

$$
y_i f(x_i) < 1
$$

For example:

```text
y =        [ 1, -1,  1]
scores =   [ 2,  0.2, 0.5]
```

Then:

```text
y * scores = [2, -0.2, 0.5]
```

Checking:

```text
[2, -0.2, 0.5] < 1
```

gives:

```text
[False, True, True]
```

So samples 2 and 3 violate the margin.

---

## Regularization

```python
alpha *= (1 - eta * lambda_param)
```

This gradually shrinks the model parameters.

Regularization prevents the model from growing unnecessarily large.

---

## Update Violating Samples

```python
alpha[violations] += eta / n_samples
```

Only samples violating

$$
y_i f(x_i) < 1
$$

receive an additional contribution.

---

# Algorithm Flow

```text
Training data X, y
        |
        v
Choose kernel
        |
        v
Build kernel matrix K
        |
        v
Initialize alpha = 0
        |
        v
For each iteration t
        |
        +--> calculate learning rate
        |
        +--> compute SVM scores
        |
        +--> calculate y * score
        |
        +--> identify margin violations
        |
        +--> shrink alpha
        |
        +--> update violating samples
        |
        +--> update bias
        |
        v
Return alpha and bias
```

---

# Linear vs RBF Kernel

| Linear Kernel | RBF Kernel |
|---|---|
| $K(x,z)=x^Tz$ | $K(x,z)=e^{-\gamma\|x-z\|^2}$ |
| Linear decision boundary | Nonlinear decision boundary |
| Faster | More computationally expensive |
| No $\gamma$ parameter | Requires $\gamma$ |
| Good for approximately linear data | Good for nonlinear patterns |

---

# Key Idea

The entire kernel SVM prediction can be summarized as:

$$
\boxed{
f(x)
=
\sum_i
\alpha_i y_i K(x_i,x)
+
b
}
$$

and classification is:

$$
\boxed{
\hat y
=
\operatorname{sign}(f(x))
}
$$

Pegasos repeatedly adjusts the model so that ideally:

$$
\boxed{
y_i f(x_i) \ge 1
}
$$

for the training samples.
