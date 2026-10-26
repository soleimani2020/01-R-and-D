# Calculate Performance Metrics for a Classification Model

> 📊 An introduction to evaluating binary classification models using confusion matrices and performance metrics.

## 📌 Overview

Classification metrics help evaluate how well a machine learning model distinguishes between positive and negative classes.

This project focuses on calculating five important performance measures:

- Confusion Matrix
- Accuracy
- F1 Score
- Specificity
- Negative Predictive Value (NPV)

The metrics are calculated from actual and predicted binary labels.

## 🎯 Problem Definition

The objective is to calculate the performance metrics using two lists of equal length:

- **Actual labels:** The ground-truth classes.
- **Predicted labels:** The classes predicted by the model.

The labels are defined as:

- `1` — Positive class
- `0` — Negative class

Both lists must contain the same number of observations.

## 🧮 Confusion Matrix

The confusion matrix summarizes the relationship between actual and predicted labels.

| Actual / Predicted | Predicted Negative (0) | Predicted Positive (1) |
|---|---:|---:|
| Actual Negative (0) | True Negative (TN) | False Positive (FP) |
| Actual Positive (1) | False Negative (FN) | True Positive (TP) |

The four outcomes are:

- **True Positive (TP):** The model correctly predicts a positive class.
- **True Negative (TN):** The model correctly predicts a negative class.
- **False Positive (FP):** The model predicts positive when the actual class is negative.
- **False Negative (FN):** The model predicts negative when the actual class is positive.

These four values are used to calculate the performance metrics.

## 📈 Performance Metrics

### 1. Accuracy

Accuracy measures the proportion of all predictions that are correct.

$$
\text{Accuracy} =
\frac{TP + TN}
{TP + TN + FP + FN}
$$

A higher accuracy indicates that a larger proportion of predictions are correct. However, accuracy can be misleading when the classes are highly imbalanced.

### 2. Precision

Precision measures the proportion of positive predictions that are actually positive.

$$
\text{Precision} =
\frac{TP}{TP + FP}
$$

High precision means that the model produces relatively few false-positive predictions.

Precision is used to calculate the F1 Score.

### 3. Recall

Recall measures the proportion of actual positive samples correctly identified by the model.

$$
\text{Recall} =
\frac{TP}{TP + FN}
$$

High recall means that the model misses relatively few actual positive samples.

Recall is particularly important when failing to identify a positive case has significant consequences.

### 4. F1 Score

The F1 Score is the harmonic mean of precision and recall. It combines both metrics into a single measure.

$$
F_1 =
2 \times
\frac{\text{Precision} \times \text{Recall}}
{\text{Precision} + \text{Recall}}
$$

An equivalent formula expressed using the confusion matrix values is:

$$
F_1 =
\frac{2TP}{2TP + FP + FN}
$$

The F1 Score is useful when both precision and recall matter, particularly when evaluating imbalanced classification problems.

### 5. Specificity

Specificity measures the proportion of actual negative samples correctly identified by the model.

$$
\text{Specificity} =
\frac{TN}{TN + FP}
$$

High specificity means that the model correctly recognizes most negative samples and produces relatively few false positives.

Specificity is also known as the **true negative rate**.

### 6. Negative Predictive Value (NPV)

Negative Predictive Value measures the proportion of negative predictions that are actually negative.

$$
\text{NPV} =
\frac{TN}{TN + FN}
$$

High NPV means that a negative prediction is likely to be correct.

Unlike specificity, which is calculated using actual negative samples, NPV is calculated using all predicted negative samples.

## 📋 Metric Summary

| Metric | Formula | What It Measures |
|---|---|---|
| Accuracy | \(\frac{TP+TN}{TP+TN+FP+FN}\) | Overall proportion of correct predictions |
| Precision | \(\frac{TP}{TP+FP}\) | Reliability of positive predictions |
| Recall | \(\frac{TP}{TP+FN}\) | Proportion of actual positives identified |
| F1 Score | \(\frac{2TP}{2TP+FP+FN}\) | Balance between precision and recall |
| Specificity | \(\frac{TN}{TN+FP}\) | Proportion of actual negatives identified |
| NPV | \(\frac{TN}{TN+FN}\) | Reliability of negative predictions |

## ⚠️ Important Considerations

- **Equal input lengths:** The actual and predicted label lists must have the same length.
- **Binary labels:** The inputs use `0` and `1` to represent negative and positive classes.
- **Zero denominators:** Some metrics are undefined when their denominators are zero. For example, precision is undefined when there are no positive predictions. The implementation must decide how to handle these cases.
- **Class imbalance:** Accuracy alone may not adequately represent model performance when one class is much more common than the other.
- **Metric selection:** The most appropriate metric depends on the classification problem and the consequences of false-positive and false-negative predictions.

## ✅ Key Takeaways

- The confusion matrix provides the foundation for evaluating binary classification.
- Accuracy measures overall correctness.
- Precision and recall focus on positive predictions and actual positive samples, respectively.
- The F1 Score combines precision and recall.
- Specificity evaluates the identification of actual negative samples.
- NPV evaluates the reliability of negative predictions.

Understanding these metrics makes it easier to evaluate classification models and choose appropriate performance measures for different machine learning applications.
