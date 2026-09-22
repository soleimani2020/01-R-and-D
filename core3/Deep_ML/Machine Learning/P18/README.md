# 📌 Lasso Regression using ISTA (from Scratch)

This project implements **Lasso Regression (L1 Regularization)** from scratch using **NumPy**, based on the **ISTA (Iterative Shrinkage-Thresholding Algorithm)**.

It demonstrates how to combine:
- Gradient Descent (for MSE loss)
- Proximal Operators (for L1 regularization)

---

## 🚀 Features

- Fully vectorized NumPy implementation  


```
# 📌 Lasso Regression using ISTA (from Scratch)

This project implements **Lasso Regression (L1 Regularization)** from scratch using **NumPy**, based on the **ISTA (Iterative Shrinkage-Thresholding Algorithm)**.

It demonstrates how to combine:

* Gradient Descent for the MSE loss
* Proximal Operators for L1 regularization
* Soft-thresholding for sparse feature selection

---

## 📘 Linear Regression

For linear regression:

```math
\hat{y} = Xw + b
```

The Mean Squared Error loss is:

```math
\mathcal{L}_{\mathrm{MSE}}(w)
=
\frac{1}{2n}
\|Xw-y\|_2^2
```

Regularization adds a penalty term to control the size of the model coefficients and reduce overfitting.

---

# 🔵 Ridge Regression — L2 Regularization

Ridge Regression adds an **L2 penalty** to the ordinary regression loss.

```math
\mathcal{L}_{\mathrm{Ridge}}(w)
=
\frac{1}{2n}
\|Xw-y\|_2^2
+
\lambda \|w\|_2^2
```

where:

```math
\|w\|_2^2
=
\sum_{j=1}^{p} w_j^2
```

Therefore:

```math
\mathcal{L}_{\mathrm{Ridge}}
=
\frac{1}{2n}
\sum_{i=1}^{n}
(y_i-\hat{y}_i)^2
+
\lambda
\sum_{j=1}^{p} w_j^2
```

Here:

* `λ` controls the strength of regularization.
* Larger `λ` means stronger penalization.
* Ridge shrinks large coefficients toward zero.
* Ridge usually does **not** make coefficients exactly zero.

Example:

```text
Before Ridge:
w = [5.2, -3.7, 0.8, 2.1]

After Ridge:
w = [3.1, -2.2, 0.4, 1.2]
```

---

## Ridge Gradient

The derivative of the L2 penalty is:

```math
\nabla_w \|w\|_2^2 = 2w
```

So the gradient becomes:

```math
\nabla_w \mathcal{L}
=
\frac{1}{n}X^T(Xw-y)
+
2\lambda w
```

Because the Ridge objective is differentiable, standard Gradient Descent can be used.

---

# 🟣 Lasso Regression — L1 Regularization

Lasso Regression adds an **L1 penalty**:

```math
\mathcal{L}_{\mathrm{Lasso}}(w)
=
\frac{1}{2n}
\|Xw-y\|_2^2
+
\lambda \|w\|_1
```

where:

```math
\|w\|_1
=
\sum_{j=1}^{p}|w_j|
```

Therefore:

```math
\mathcal{L}_{\mathrm{Lasso}}
=
\frac{1}{2n}
\sum_{i=1}^{n}
(y_i-\hat{y}_i)^2
+
\lambda
\sum_{j=1}^{p}|w_j|
```

Unlike Ridge, Lasso can force some coefficients to become exactly zero.

Example:

```text
Before Lasso:
w = [5.2, -3.7, 0.8, 2.1]

After Lasso:
w = [3.8, -2.5, 0.0, 0.0]
```

This means Lasso performs a form of **automatic feature selection**.

---

# ⚠️ Why Lasso Is Different

The L1 penalty contains the absolute value:

```math
|w|
```

Its derivative is:

```math
\frac{d|w|}{dw}
=
\begin{cases}
1, & w>0 \\
-1, & w<0
\end{cases}
```

At:

```math
w=0
```

the ordinary derivative is not defined.

Because of this non-smooth point, Lasso is commonly optimized using methods such as:

* ISTA
* FISTA
* Coordinate Descent

---

# 🔁 ISTA

ISTA stands for:

**Iterative Shrinkage-Thresholding Algorithm**

It separates the Lasso objective into two parts:

```math
f(w)
=
\frac{1}{2n}
\|Xw-y\|_2^2
+
\lambda\|w\|_1
```

The first term is smooth:

```math
\frac{1}{2n}
\|Xw-y\|_2^2
```

The second term is non-smooth:

```math
\lambda\|w\|_1
```

ISTA handles these two parts separately.

---

## Step 1 — Gradient Descent

First, take a Gradient Descent step using only the MSE part:

```math
z
=
w^{(k)}
-
\eta
\nabla\mathcal{L}_{\mathrm{MSE}}
\left(w^{(k)}\right)
```

where:

```math
\nabla\mathcal{L}_{\mathrm{MSE}}
=
\frac{1}{n}
X^T(Xw-y)
```

Here:

* `η` is the learning rate
* `z` is the intermediate parameter vector

---

## Step 2 — Soft Thresholding

After the gradient step, apply the soft-thresholding operator:

```math
w^{(k+1)}
=
S_{\eta\lambda}(z)
```

The soft-thresholding operator is:

```math
S_{\tau}(z)
=
\begin{cases}
z-\tau, & z>\tau \\
0, & |z|\le\tau \\
z+\tau, & z<-\tau
\end{cases}
```

where:

```math
\tau = \eta\lambda
```

This step shrinks the coefficients toward zero.

If a coefficient is already small enough, it becomes exactly zero.

---

## Example of Soft Thresholding

Suppose:

```text
τ = 0.2
```

Then:

```text
1.5  →  1.3
-0.8 → -0.6
0.1  →  0.0
```

So the operator performs both:

* coefficient shrinkage
* automatic feature selection

---

# 🆚 Ridge vs. Lasso

| Property             | Ridge            | Lasso                     |
| -------------------- | ---------------- | ------------------------- |
| Regularization       | L2               | L1                        |
| Penalty              | `λ Σ w²`         | `λ Σ \|w\|`               |
| Shrinks coefficients | ✅ Yes            | ✅ Yes                     |
| Produces exact zeros | Usually no       | ✅ Yes                     |
| Feature selection    | ❌ No             | ✅ Yes                     |
| Objective smooth     | ✅ Yes            | ❌ Not at zero             |
| Common optimization  | Gradient Descent | ISTA / Coordinate Descent |

---

# 🧠 Key Idea

A simple way to remember the difference:

```text
Ridge
↓
Shrink coefficients
```

```text
Lasso
↓
Shrink coefficients
+
Set some coefficients exactly to zero
```

For ISTA:

```text
Gradient Descent
        +
Soft Thresholding
        =
Lasso Optimization
```

---

## 🚀 Features

* Fully vectorized NumPy implementation
* Lasso Regression from scratch
* ISTA optimization
* Soft-thresholding operator
* L1 regularization
* Automatic feature selection
* No machine-learning library required

---

## 📦 Requirements

```bash
pip install numpy
```

---

## 🎯 Main Concepts

This project demonstrates:

* Linear Regression
* Mean Squared Error
* Ridge Regression
* Lasso Regression
* L1 and L2 Regularization
* Gradient Descent
* Proximal Gradient Methods
* Soft Thresholding
* Sparse Models
* Feature Selection
