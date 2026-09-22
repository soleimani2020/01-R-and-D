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

This project implements **Lasso Regression** from scratch using **NumPy** and the **ISTA (Iterative Shrinkage-Thresholding Algorithm)**.

It combines:

* Gradient Descent
* L1 Regularization
* Soft Thresholding

---

## 📘 Linear Regression

A linear model predicts:

```math
\hat{y} = Xw
```

where:

* `X` = input features
* `w` = model weights
* `ŷ` = predicted values

The Mean Squared Error is:

```math
L = \frac{1}{2n}\sum_{i=1}^{n}(y_i-\hat{y}_i)^2
```

---

# 🔵 Ridge Regression

Ridge adds an **L2 penalty**:

```math
L =
\frac{1}{2n}\sum_{i=1}^{n}(y_i-\hat{y}_i)^2
+
\lambda\sum_{j=1}^{p}w_j^2
```

The penalty is:

```math
\lambda\sum w_j^2
```

This makes large weights expensive.

So Ridge pushes the weights toward zero:

```text
5.0   →   3.2
-4.0  →  -2.5
1.0   →   0.6
```

But Ridge usually does not make weights exactly zero.

### Ridge Gradient

The gradient is:

```math
\nabla L =
\frac{1}{n}X^T(Xw-y)
+
2\lambda w
```

So Ridge can be optimized with ordinary Gradient Descent.

---

# 🟣 Lasso Regression

Lasso adds an **L1 penalty**:

```math
L =
\frac{1}{2n}\sum_{i=1}^{n}(y_i-\hat{y}_i)^2
+
\lambda\sum_{j=1}^{p}|w_j|
```

The penalty is:

```math
\lambda\sum |w_j|
```

Lasso also shrinks the weights, but it can make some of them exactly zero.

Example:

```text
Before:
w = [4.0, -2.0, 0.3, 0.1]

After:
w = [3.2, -1.2, 0.0, 0.0]
```

This means Lasso can automatically remove unimportant features.

---

# ⚠️ Why Lasso Needs ISTA

The L1 penalty contains:

```math
|w|
```

This function is not differentiable at:

```math
w = 0
```

So instead of applying ordinary Gradient Descent to the entire loss, ISTA uses two steps.

---

# 🔁 ISTA

## Step 1 — Gradient Step

First, update the weights using only the MSE part:

```math
z = w - \eta \nabla L
```

where:

* `η` = learning rate
* `z` = temporary updated weight

The MSE gradient is:

```math
\nabla L =
\frac{1}{n}X^T(Xw-y)
```

---

## Step 2 — Soft Thresholding

After the gradient step, shrink each value of `z`.

```math
w_{\text{new}} = S(z)
```

A simple form of the soft-thresholding rule is:

```math
S(z) =
\begin{cases}
z-\tau, & z>\tau \\
0, & -\tau \le z \le \tau \\
z+\tau, & z<-\tau
\end{cases}
```

where:

```math
\tau = \eta\lambda
```

So:

```text
if z is large positive:
    subtract τ

if z is large negative:
    add τ

if z is close to zero:
    set it to zero
```

Example:

```text
τ = 0.2

1.5   →   1.3
-0.8  →  -0.6
0.1   →   0.0
```

This is the key idea behind Lasso.

---

# 🆚 Ridge vs. Lasso

| Property              | Ridge            | Lasso |   |   |
| --------------------- | ---------------- | ----- | - | - |
| Penalty               | `w²`             | `     | w | ` |
| Regularization        | L2               | L1    |   |   |
| Shrinks weights       | Yes              | Yes   |   |   |
| Can make weights zero | Usually no       | Yes   |   |   |
| Feature selection     | No               | Yes   |   |   |
| Common method         | Gradient Descent | ISTA  |   |   |

---

# 🧠 Easy Way to Remember

```text
Ridge
→ shrink weights
```

```text
Lasso
→ shrink weights
→ set small weights to zero
```

And ISTA is simply:

```text
Gradient Step
      +
Soft Thresholding
      =
Lasso Optimization
```
