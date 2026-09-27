# Binary Classification with Logistic Regression

A concise guide to logistic regression for binary classification.

---

## Table of Contents

- [Overview](#overview)
- [Model](#model)
- [Sigmoid Function](#sigmoid-function)
- [Classification Threshold](#classification-threshold)
- [Decision Boundary](#decision-boundary)
- [Log-Odds Interpretation](#log-odds-interpretation)
- [Training Objective](#training-objective)
- [Python Example](#python-example)
- [Summary](#summary)

---

## Overview

Logistic Regression is a linear model used for **binary classification**, where the target has two possible classes:

$$
y \in \{0, 1\}
$$

The model computes a linear score:

$$
z = Xw + b
$$

where:

- $X$ = feature matrix
- $w$ = model weights
- $b$ = bias / intercept
- $z$ = linear score

The linear score $z$ can be any real number, so it is not yet a probability.

---

## Sigmoid Function

To convert the linear score $z$ into a probability between `0` and `1`, we apply the **sigmoid function**:

$$
\sigma(z) = \frac{1}{1 + e^{-z}}
$$

Therefore, the probability that the class is `1` is:

$$
P(y=1 \mid X) = \sigma(Xw + b)
$$

The probability that the class is `0` is:

$$
P(y=0 \mid X) = 1 - P(y=1 \mid X)
$$

### Key Properties

$$
0 < \sigma(z) < 1
$$

$$
\sigma(0) = 0.5
$$

$$
\sigma(z) \to 1 \quad \text{as } z \to \infty
$$

$$
\sigma(z) \to 0 \quad \text{as } z \to -\infty
$$

The sigmoid function is monotonic, meaning larger $z$ values produce larger probabilities.

---

## Classification Threshold

After computing the probability, we convert it into a class prediction.

Using the default threshold of `0.5`:

$$
\hat{y}
=
\begin{cases}
1, & \text{if } P(y=1 \mid X) \ge 0.5 \\
0, & \text{if } P(y=1 \mid X) < 0.5
\end{cases}
$$

Equivalently:

```text
Probability >= 0.5  -> Class 1
Probability <  0.5  -> Class 0
