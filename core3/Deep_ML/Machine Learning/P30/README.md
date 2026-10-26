# Binary Classification with Logistic Regression

> 📊 **Machine Learning** · 🎯 **Binary Classification** · 📈 **Probability Estimation** · 🧮 **Linear Models**

## 📌 Project Overview

**Logistic Regression** is a supervised machine learning algorithm used for binary classification. It estimates the probability that an observation belongs to one of two classes and converts that probability into a class prediction using a decision threshold.

Although its name includes *regression*, logistic regression is commonly used for classification tasks.

This guide covers the model's mathematical foundations, the sigmoid function, decision boundaries, log-odds interpretation, and training objective.

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

## 1. Overview

Logistic regression is designed for classification problems in which the target variable has two possible values:

\[
y\in\{0,1\}
\]

The model first calculates a linear score:

\[
z=Xw+b
\]

where:

- \(X\) represents the input features.
- \(w\) represents the model weights.
- \(b\) is the bias or intercept.
- \(z\) is the linear score.

The linear score can take any real value, so it cannot directly be interpreted as a probability. The sigmoid function transforms it into a value between 0 and 1.

## 2. The Model

For a single observation \(x\), logistic regression calculates:

\[
z=x^\top w+b
\]

For a dataset containing multiple observations, the corresponding vectorized expression is:

\[
z=Xw+b
\]

The resulting score is passed through the sigmoid function to estimate the probability of the positive class.

## 3. The Sigmoid Function

The sigmoid function converts a real-valued score into a probability:

\[
\boxed{
\sigma(z)=\frac{1}{1+e^{-z}}
}
\]

The probability that an observation belongs to class 1 is:

\[
P(y=1\mid x)=\sigma(x^\top w+b)
\]

Writing this probability as \(p\), the probability of class 0 is:

\[
P(y=0\mid x)=1-p
\]

### Key Properties

The sigmoid function has several important properties:

\[
0<\sigma(z)<1
\]

\[
\sigma(0)=0.5
\]

\[
\lim_{z\to+\infty}\sigma(z)=1
\]

\[
\lim_{z\to-\infty}\sigma(z)=0
\]

The function is monotonically increasing: larger linear scores produce larger estimated probabilities.

## 4. Classification Threshold

Logistic regression estimates probabilities, but a classification task usually requires a discrete prediction.

Using the default threshold of 0.5:

\[
\hat y=
\begin{cases}
1, & p\geq 0.5\\
0, & p<0.5
\end{cases}
\]

In practice:

```text
Probability >= 0.5  -> Class 1
Probability <  0.5  -> Class 0
```

The threshold can be adjusted depending on the application. For example, a lower threshold may increase recall at the expense of precision.

The threshold of 0.5 is a common default, not a universally optimal choice.

## 5. Decision Boundary

Because the sigmoid function is monotonically increasing and \(\sigma(0)=0.5\), the default classification rule can also be expressed using the linear score:

\[
\hat y=
\begin{cases}
1, & x^\top w+b\geq 0\\
0, & x^\top w+b<0
\end{cases}
\]

The decision boundary is defined by:

\[
\boxed{x^\top w+b=0}
\]

This boundary separates the feature space into regions associated with the two predicted classes.

For two input features, the boundary is generally a straight line. With more features, it becomes a hyperplane.

Although the estimated probabilities are produced by a nonlinear sigmoid transformation, the standard logistic regression decision boundary is linear in the original feature space.

## 6. Log-Odds Interpretation

Logistic regression can also be understood through the concepts of odds and log-odds.

The odds of class 1 are:

\[
\text{Odds}=\frac{p}{1-p}
\]

Taking the natural logarithm gives the log-odds:

\[
\log\left(\frac{p}{1-p}\right)=x^\top w+b
\]

Therefore, logistic regression models the **log-odds of the positive class as a linear function of the input features**.

### Interpreting the Coefficients

For a feature \(x_j\), its coefficient \(w_j\) represents the change in log-odds associated with a one-unit increase in that feature, holding the other features constant.

Exponentiating the coefficient gives:

\[
e^{w_j}
\]

This is the corresponding multiplicative change in the odds for a one-unit increase in the feature, assuming the model's other inputs remain fixed.

A positive coefficient increases the estimated log-odds as the feature increases, while a negative coefficient decreases them.

## 7. Training Objective

The model parameters \(w\) and \(b\) are learned by minimizing the **binary cross-entropy loss**, also known as log loss.

For \(N\) training observations:

\[
\boxed{
J(w,b)=
-\frac{1}{N}\sum_{i=1}^{N}
\left[
y_i\log(p_i)
+(1-y_i)\log(1-p_i)
\right]
}
\]

where:

\[
p_i=\sigma(x_i^\top w+b)
\]

The loss measures how well the predicted probabilities match the true labels.

- Confident, correct predictions receive a low loss.
- Uncertain predictions generally receive a higher loss.
- Confident, incorrect predictions receive a particularly high loss.

Minimizing binary cross-entropy encourages the model to assign high probabilities to actual positive observations and low probabilities to actual negative observations.

Regularization can also be added to the objective to reduce overfitting.

## 8. Python Example

The following example demonstrates logistic regression using `scikit-learn`.

```python
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, classification_report

# Generate a binary classification dataset
X, y = make_classification(
    n_samples=1000,
    n_features=10,
    n_informative=6,
    n_redundant=2,
    random_state=42
)

# Split the data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.2,
    random_state=42,
    stratify=y
)

# Scale the features and train the model
model = make_pipeline(
    StandardScaler(),
    LogisticRegression(max_iter=1000)
)

model.fit(X_train, y_train)

# Predict class labels and probabilities
y_pred = model.predict(X_test)
y_proba = model.predict_proba(X_test)[:, 1]

# Evaluate the model
print("Accuracy:", accuracy_score(y_test, y_pred))
print(classification_report(y_test, y_pred))
```

### Understanding the Example

1. **Data generation:** creates a synthetic binary classification dataset.
2. **Train-test split:** reserves 20% of the observations for testing.
3. **Feature scaling:** standardizes the features using statistics learned from the training data.
4. **Model training:** fits logistic regression to the training data.
5. **Class prediction:** produces discrete class labels.
6. **Probability estimation:** returns the estimated probability of class 1.
7. **Evaluation:** measures performance on the held-out test set.

The example is illustrative; actual performance depends on the dataset and its underlying patterns.

## 9. Summary

The logistic regression workflow can be summarized in four steps.

**Step 1 — Calculate the linear score:**

\[
z=x^\top w+b
\]

**Step 2 — Convert the score into a probability:**

\[
p=\sigma(z)=\frac{1}{1+e^{-z}}
\]

**Step 3 — Apply the classification threshold:**

```text
p >= 0.5  -> Predict class 1
p <  0.5  -> Predict class 0
```

**Step 4 — Learn the parameters:**

Minimize binary cross-entropy to obtain weights and a bias that fit the training data.

### 📌 Key Takeaway

Logistic regression is an interpretable linear classifier that estimates class probabilities using the sigmoid function. Its parameters are learned through an optimization objective, and its predictions are determined by applying a classification threshold to the estimated probabilities.
