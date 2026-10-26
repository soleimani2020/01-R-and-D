# Binary Classification with Logistic Regression

> 📊 **Machine Learning** · 🎯 **Binary Classification** · 📈 **Probability Estimation** · 🧮 **Linear Models**

## 📌 Project Overview

**Logistic Regression** is a supervised machine learning algorithm used for binary classification. It estimates the probability that an observation belongs to one of two classes and converts that probability into a class prediction using a decision threshold.

Although its name includes *regression*, logistic regression is commonly used for classification tasks.

This guide covers the mathematical foundations of the model, the sigmoid function, classification thresholds, decision boundaries, log-odds interpretation, and the training objective.

## 📑 Table of Contents

1. [Overview](#1-overview)
2. [The Model](#2-the-model)
3. [The Sigmoid Function](#3-the-sigmoid-function)
4. [Classification Threshold](#4-classification-threshold)
5. [Decision Boundary](#5-decision-boundary)
6. [Log-Odds Interpretation](#6-log-odds-interpretation)
7. [Training Objective](#7-training-objective)
8. [Python Example](#8-python-example)
9. [Summary](#9-summary)

---

## 1. Overview

Logistic regression is designed for binary classification problems in which the target variable has two possible values:

$$
y \in \{0,1\}
$$

The model first calculates a linear score:

$$
z = Xw + b
$$

where:

- `X` represents the input feature matrix.
- `w` represents the model weights.
- `b` is the bias or intercept.
- `z` is the linear score.

The linear score can take any real value, so it cannot directly be interpreted as a probability. The sigmoid function transforms it into a value between 0 and 1.

---

## 2. The Model

For a single observation, the model calculates:

$$
z = x^\top w + b
$$

For an entire dataset, the vectorized expression is:

$$
z = Xw + b
$$

The resulting score is passed through the sigmoid function to estimate the probability of the positive class.

The model learns the weights and bias from the training data by minimizing a suitable loss function.

---

## 3. The Sigmoid Function

The sigmoid function converts a real-valued score into a probability:

$$
\sigma(z) = \frac{1}{1 + e^{-z}}
$$

The probability that an observation belongs to class 1 is:

$$
P(y = 1 \mid x) = \sigma(x^\top w + b)
$$

If we denote this probability by `p`, the probability of class 0 is:

$$
P(y = 0 \mid x) = 1 - p
$$

### Key Properties

The sigmoid function has the following properties:

**Output range:**

$$
0 < \sigma(z) < 1
$$

**At zero:**

$$
\sigma(0) = 0.5
$$

**As the score approaches positive infinity:**

$$
\lim_{z \to +\infty} \sigma(z) = 1
$$

**As the score approaches negative infinity:**

$$
\lim_{z \to -\infty} \sigma(z) = 0
$$

The sigmoid function is monotonically increasing. Therefore, a larger linear score produces a larger estimated probability of class 1.

---

## 4. Classification Threshold

Logistic regression estimates probabilities, but classification generally requires a discrete prediction.

Using the default threshold of `0.5`:

$$
\hat{y} =
\begin{cases}
1 & \text{if } p \geq 0.5 \\
0 & \text{if } p < 0.5
\end{cases}
$$

where:

$$
p = P(y = 1 \mid x)
$$

The rule can also be written as:

```text
Probability >= 0.5  -> Class 1
Probability <  0.5  -> Class 0
```

The threshold can be adjusted according to the application.

For example, lowering the threshold may increase recall while reducing precision. The appropriate threshold depends on the relative costs of false positives and false negatives.

---

## 5. Decision Boundary

Because the sigmoid function is monotonically increasing and:

$$
\sigma(0) = 0.5
$$

the classification rule can also be expressed directly in terms of the linear score.

$$
\hat{y} =
\begin{cases}
1 & \text{if } x^\top w + b \geq 0 \\
0 & \text{if } x^\top w + b < 0
\end{cases}
$$

The decision boundary is defined by:

$$
\boxed{x^\top w + b = 0}
$$

This boundary separates the feature space into regions associated with the two predicted classes.

- With one feature, the boundary is a point on the number line.
- With two features, the boundary is generally a straight line.
- With more features, the boundary is a hyperplane.

Although the sigmoid function is nonlinear, standard logistic regression produces a linear decision boundary in the original feature space.

---

## 6. Log-Odds Interpretation

Logistic regression can also be understood through odds and log-odds.

### Odds

The odds of class 1 are defined as:

$$
\text{Odds} = \frac{p}{1-p}
$$

where `p` is the probability of class 1.

### Log-Odds

Taking the natural logarithm gives:

$$
\log\left(\frac{p}{1-p}\right) = x^\top w + b
$$

Therefore, logistic regression models the **log-odds of the positive class as a linear function of the input features**.

### Interpreting the Coefficients

For a feature `x_j`, the coefficient `w_j` represents the change in log-odds associated with a one-unit increase in that feature, holding the other features constant.

Exponentiating the coefficient gives the odds ratio:

$$
\text{Odds Ratio} = e^{w_j}
$$

Interpretation:

- If `w_j > 0`, increasing the feature increases the log-odds of class 1.
- If `w_j < 0`, increasing the feature decreases the log-odds of class 1.
- If `w_j = 0`, the feature has no direct contribution to the linear score.

The odds ratio describes a multiplicative change in odds, not a constant change in probability.

---

## 7. Training Objective

The model parameters `w` and `b` are learned by minimizing **binary cross-entropy loss**, also known as log loss.

For `N` training observations:

$$
J(w,b) =
-\frac{1}{N}
\sum_{i=1}^{N}
\left[
y_i \log(p_i)
+
(1-y_i)\log(1-p_i)
\right]
$$

where the predicted probability for observation `i` is:

$$
p_i = \frac{1}{1 + e^{-(x_i^\top w+b)}}
$$

The loss measures how well the predicted probabilities match the true labels.

- Confident, correct predictions receive a low loss.
- Uncertain predictions generally receive a higher loss.
- Confident, incorrect predictions receive a particularly high loss.

Minimizing binary cross-entropy encourages the model to assign high probabilities to positive observations and low probabilities to negative observations.

Regularization can also be added to the objective to reduce overfitting.

---

## 8. Python Example

The following example demonstrates logistic regression using `scikit-learn`.

### Import Libraries

```python
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, classification_report
```

### Generate a Dataset

```python
X, y = make_classification(
    n_samples=1000,
    n_features=10,
    n_informative=6,
    n_redundant=2,
    random_state=42
)
```

This creates a synthetic binary classification dataset.

### Split the Data

```python
X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.2,
    random_state=42,
    stratify=y
)
```

The dataset is divided into training and test sets. The `stratify=y` argument helps preserve the class proportions in both sets.

### Train the Model

```python
model = make_pipeline(
    StandardScaler(),
    LogisticRegression(max_iter=1000)
)

model.fit(X_train, y_train)
```

The pipeline standardizes the features using training-set statistics before fitting the logistic regression model.

### Make Predictions

```python
y_pred = model.predict(X_test)

y_proba = model.predict_proba(X_test)[:, 1]
```

- `y_pred` contains predicted class labels.
- `y_proba` contains the estimated probabilities of class 1.

### Evaluate the Model

```python
print("Accuracy:", accuracy_score(y_test, y_pred))
print(classification_report(y_test, y_pred))
```

These metrics summarize the model's classification performance on the held-out test set.

The example is illustrative; actual performance depends on the dataset and its underlying patterns.

---

## 9. Summary

Logistic regression follows four main steps.

### Step 1 — Calculate the Linear Score

$$
z = x^\top w + b
$$

### Step 2 — Convert the Score into a Probability

$$
p = \sigma(z) = \frac{1}{1 + e^{-z}}
$$

### Step 3 — Apply the Classification Threshold

$$
\hat{y} =
\begin{cases}
1 & \text{if } p
