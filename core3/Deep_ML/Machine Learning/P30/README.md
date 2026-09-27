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

## Model

The model first calculates a linear score:

$$
z = Xw + b
$$

Then it converts this score into a probability using the sigmoid function.

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
```

---

## Decision Boundary

Because the sigmoid function is monotonic and $\sigma(0) = 0.5$, the threshold rule can also be written directly in terms of the linear score:

$$
\hat{y}
=
\begin{cases}
1, & \text{if } Xw + b \ge 0 \\
0, & \text{if } Xw + b < 0
\end{cases}
$$

The boundary where the model is exactly undecided is:

$$
Xw + b = 0
$$

This is called the **decision boundary**.

---

## Log-Odds Interpretation

Logistic regression can also be interpreted using odds.

The odds of class `1` are:

$$
\frac{P(y=1 \mid X)}{P(y=0 \mid X)}
=
\frac{p}{1-p}
$$

Taking the log gives:

$$
\log\left(\frac{p}{1-p}\right) = Xw + b
$$

So logistic regression models the **log-odds** of the positive class as a linear function of the input features.

---

## Training Objective

The parameters $w$ and $b$ are learned by minimizing the **binary cross-entropy loss**, also called the log loss:

$$
J(w,b)
=
-\frac{1}{N}
\sum_{i=1}^{N}
\left[
y_i \log(p_i) + (1-y_i)\log(1-p_i)
\right]
$$

where:

$$
p_i = \sigma(x_i w + b)
$$

This loss penalizes incorrect confident predictions more strongly than uncertain incorrect predictions.

---


## Summary

1. Compute the linear score:

$$
z = Xw + b
$$

2. Convert it into a probability using the sigmoid function:

$$
p = \sigma(z)
$$

3. Apply a threshold, usually `0.5`:

```text
p >= 0.5  -> predict class 1
p <  0.5  -> predict class 0
```

4. Equivalently, because sigmoid is monotonic:

```text
z >= 0  -> predict class 1
z <  0  -> predict class 0
```

Logistic regression is therefore a simple, interpretable linear classifier that outputs probabilities and makes binary decisions through a threshold.
