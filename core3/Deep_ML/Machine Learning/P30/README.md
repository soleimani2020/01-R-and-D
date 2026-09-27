# Binary Classification with Logistic Regression

## Overview

Logistic Regression is used for **binary classification**, where the target has two possible classes:

$$
y \in \{0,1\}
$$

The model first calculates a linear score:

$$
z = Xw + b
$$

where:

- $X$ = feature matrix
- $w$ = model weights
- $b$ = bias/intercept
- $z$ = linear score

---

## 1. Sigmoid Function

The linear score $z$ can be any real number.

To convert it into a probability between `0` and `1`, we use the sigmoid function:

$$
\sigma(z)
=
\frac{1}{1 + e^{-z}}
$$

Therefore:

$$
P(y=1 \mid X)
=
\sigma(Xw+b)
$$

The output satisfies:

$$
0 < \sigma(z) < 1
$$

---

## 2. Classification Threshold

After calculating the probability, we convert it into a class prediction.

Using a threshold of `0.5`:

$$
\hat{y}
=
\begin{cases}
1, & \text{if } P(y=1 \mid X) \ge 0.5 \\
0, & \text{if } P(y=1 \mid X) < 0.5
\end{cases}
$$

So:

```text
Probability >= 0.5  -> Class 1
Probability <  0.5  -> Class 0
