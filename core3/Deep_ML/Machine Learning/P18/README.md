# 📌 Lasso Regression using ISTA (from Scratch)

This project implements **Lasso Regression (L1 Regularization)** from scratch using **NumPy**, based on the **ISTA (Iterative Shrinkage-Thresholding Algorithm)**.

It demonstrates how to combine:
- Gradient Descent (for MSE loss)
- Proximal Operators (for L1 regularization)

---

## 🚀 Features

- Fully vectorized NumPy implementation  

## 📘 Ridge and Lasso Regression

Both **Ridge** and **Lasso** are regularized versions of linear regression.

For ordinary linear regression:

```math
\hat{y} = Xw + b
```

and the Mean Squared Error loss is:

```math
\mathcal{L}_{\text{MSE}}(w)
=
\frac{1}{2n}
\|Xw-y\|_2^2
```

Regularization adds a penalty term to reduce overfitting and control the size of the coefficients.

---

## 🔵 Ridge Regression — L2 Regularization

Ridge adds the squared magnitude of the coefficients:

```math
\mathcal{L}_{\text{Ridge}}(w)
=
\frac{1}{2n}\|Xw-y\|_2^2
+
\lambda \|w\|_2^2
```

where:

```math
\|w\|_2^2
=
\sum_{j=1}^{p} w_j^2
```

So explicitly:

```math
\mathcal{L}_{\text{Ridge}}
=
\frac{1}{2n}
\sum_{i=1}^{n}
(y_i-\hat{y}_i)^2
+
\lambda
\sum_{j=1}^{p} w_j^2
```

The parameter `λ` controls the strength of the regularization.

* If `λ = 0`, Ridge becomes ordinary linear regression.
* If `λ` increases, large coefficients are penalized more strongly.
* Ridge usually makes coefficients smaller, but not exactly zero.

Example:

```text
Before Ridge:
w = [5.2, -3.7, 0.8, 2.1]

After Ridge:
w = [3.1, -2.2, 0.4, 1.2]
```

### Ridge Gradient

The derivative of the L2 penalty is:

```math
\nabla_w \|w\|_2^2 = 2w
```

Therefore:

```math
\nabla_w \mathcal{L}
=
\frac{1}{n}X^T(Xw-y)
+
2\lambda w
```

Because Ridge is differentiable, it can be optimized directly using gradient descent.

---

## 🟣 Lasso Regression — L1 Regularization

Lasso adds the absolute values of the coefficients:

```math
\mathcal{L}_{\text{Lasso}}(w)
=
\frac{1}{2n}\|Xw-y\|_2^2
+
\lambda \|w\|_1
```

where:

```math
\|w\|_1
=
\sum_{j=1}^{p}|w_j|
```

So explicitly:

```math
\mathcal{L}_{\text{Lasso}}
=
\frac{1}{2n}
\sum_{i=1}^{n}
(y_i-\hat{y}_i)^2
+
\lambda
\sum_{j=1}^{p}|w_j|
```

Unlike Ridge, Lasso can force coefficients to become exactly zero.

Example:

```text
Before Lasso:
w = [5.2, -3.7, 0.8, 2.1]

After Lasso:
w = [3.8, -2.5, 0.0, 0.0]
```

This means Lasso can perform automatic feature selection.

---

## ⚠️ Why Lasso Needs a Different Optimization Method

The L1 penalty contains:

```math
|w|
```

Its derivative is:

```math
\frac{d|w|}{dw}
=
\begin{cases}
1 & w > 0 \\
-1 & w < 0
\end{cases}
```

but at:

```math
w = 0
```

the ordinary derivative is not defined.

For this reason, algorithms such as **ISTA** are commonly used.

---

## 🔁 ISTA

ISTA separates the Lasso objective into two parts:

```math
f(w)
=
\underbrace{
\frac{1}{2n}\|Xw-y\|_2^2
}_{\text{smooth MSE}}
+
\underbrace{
\lambda\|w\|_1
}_{\text{non-smooth L1 penalty}}
```

### Step 1: Gradient Descent

First, perform a gradient step on the MSE part:

```math
z
=
w^{(k)}
-
\eta
\nabla \mathcal{L}_{\text{MSE}}
\left(w^{(k)}\right)
```

with:

```math
\nabla \mathcal{L}_{\text{MSE}}
=
\frac{1}{n}
X^T(Xw-y)
```


### Step 2: Soft Thresholding

Then apply the soft-thresholding operator:

```math
w^{(k+1)}
=
S_{\eta\lambda}(z)
```

where:

```math
S_{\tau}(z) =
\begin{cases}
z-\tau, & z>\tau \\
0, & |z|\le\tau \\
z+\tau, & z<-\tau
\end{cases}




This operation shrinks large coefficients and sets sufficiently small coefficients exactly to zero.

Example with:

```text
τ = 0.2
```

```text
1.5  →  1.3
-0.8 → -0.6
0.1  →  0.0
```

---

## 🆚 Ridge vs. Lasso

| Property                         | Ridge            | Lasso                     |   |   |
| -------------------------------- | ---------------- | ------------------------- | - | - |
| Regularization                   | L2               | L1                        |   |   |
| Penalty                          | `λ Σ w²`         | `λ Σ                      | w | ` |
| Shrinks coefficients             | ✅ Yes            | ✅ Yes                     |   |   |
| Coefficients become exactly zero | Usually no       | ✅ Yes                     |   |   |
| Feature selection                | ❌ No             | ✅ Yes                     |   |   |
| Objective smooth                 | ✅ Yes            | ❌ No at zero              |   |   |
| Common optimization              | Gradient Descent | ISTA / Coordinate Descent |   |   |

A simple way to remember the difference:

```text
Ridge  → Shrink coefficients
Lasso  → Shrink coefficients + Select features
```

In this project, ISTA combines:

```text
Gradient Descent
        +
Soft Thresholding
        =
Lasso Optimization
```
