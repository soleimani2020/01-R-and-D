# Calculate Jaccard Index for Binary Classification

> 📊 **Machine Learning** · 🎯 **Binary Classification** · 📈 **Model Evaluation** · 🧮 **Classification Metrics**

## 📌 Project Overview

The **Jaccard Index**, also known as **Intersection over Union (IoU)**, is a metric used to measure the similarity between two sets.

In binary classification, it compares the overlap between the actual positive labels and the predicted positive labels.

The Jaccard Index ranges from 0 to 1:

\[
0 \leq J \leq 1
\]

- **0:** No overlap between the actual and predicted positive sets.
- **1:** Perfect agreement between the actual and predicted positive sets.

The Jaccard Index is particularly useful when evaluating classification models on imbalanced datasets because it focuses on positive predictions rather than counting true negatives.

## 🎯 Main Goal

Understand, calculate, and interpret the Jaccard Index as a metric for evaluating binary classification models.

## 🧮 1. Mathematical Formula

The general Jaccard Index is defined as:

\[
J(A,B)=\frac{|A\cap B|}{|A\cup B|}
\]

where:

- \(A\) is the set of actual positive samples.
- \(B\) is the set of predicted positive samples.
- \(A\cap B\) represents samples that are positive in both sets.
- \(A\cup B\) represents samples that are positive in at least one set.

For binary classification, the formula becomes:

\[
\boxed{
J=\frac{TP}{TP+FP+FN}
}
\]

where:

| Term | Meaning |
|---|---|
| \(TP\) | True Positives: actual positives correctly predicted as positive |
| \(FP\) | False Positives: actual negatives incorrectly predicted as positive |
| \(FN\) | False Negatives: actual positives incorrectly predicted as negative |

**Important:** True Negatives (\(TN\)) are not included in the formula.

## 📊 2. Understanding the Metric

The Jaccard Index evaluates how well the model's positive predictions match the actual positive labels.

- A high score indicates substantial overlap between actual and predicted positives.
- A low score indicates limited overlap.
- A score of 1 indicates perfect agreement on the positive class.
- A score of 0 indicates no shared positive predictions, provided the union is nonempty.

Unlike accuracy, the Jaccard Index does not increase simply because the model correctly predicts many negative samples.

## 🔢 3. Worked Example

Suppose a binary classification model produces the following results:

- True Positives: \(TP=40\)
- False Positives: \(FP=10\)
- False Negatives: \(FN=20\)

Substituting these values into the formula gives:

\[
J=\frac{40}{40+10+20}
\]

\[
J=\frac{40}{70}\approx 0.5714
\]

Therefore, the Jaccard Index is approximately **0.57**.

This indicates that the intersection of the actual and predicted positive sets represents about 57.14% of their union.

## 💻 4. Python Implementation

### Calculate the Jaccard Index manually

```python
def jaccard_index(y_true, y_pred):
    tp = sum((yt == 1 and yp == 1) for yt, yp in zip(y_true, y_pred))
    fp = sum((yt == 0 and yp == 1) for yt, yp in zip(y_true, y_pred))
    fn = sum((yt == 1 and yp == 0) for yt, yp in zip(y_true, y_pred))

    denominator = tp + fp + fn

    if denominator == 0:
        return 0.0

    return tp / denominator
```

**Note:** This implementation assumes labels are encoded as `0` and `1`, the two input arrays have equal lengths, and an empty union should return `0.0`. Other conventions may handle the empty-union case differently.

### Using scikit-learn

The Jaccard Index can also be calculated using `scikit-learn`:

```python
from sklearn.metrics import jaccard_score

score = jaccard_score(y_true, y_pred, pos_label=1)
print(f"Jaccard Index: {score:.4f}")
```

For binary classification, `pos_label=1` specifies that label `1` is the positive class.

## 🔗 5. Relationship to Other Metrics

The Jaccard Index is closely related to the F1-score.

For binary classification:

\[
F_1=\frac{2TP}{2TP+FP+FN}
\]

The relationship between the two metrics is:

\[
\boxed{
F_1=\frac{2J}{1+J}
}
\]

Both metrics focus on the positive class, but the F1-score generally has a higher numerical value for the same nonzero overlap.

## 🛠️ Technologies & Tools

| Technology | Purpose |
|---|---|
| Python | Implementation |
| NumPy | Numerical operations, if used |
| scikit-learn | Metric calculation |
| Jupyter Notebook | Interactive experimentation, if used |

## 🚀 Expected Outcomes

By completing this project, you should be able to:

- Explain the mathematical definition of the Jaccard Index.
- Calculate the metric using true positives, false positives, and false negatives.
- Implement the calculation in Python.
- Use scikit-learn to evaluate binary classification predictions.
- Interpret the result and understand why true negatives are excluded.
- Distinguish the Jaccard Index from accuracy and the F1-score.

## 📌 Key Takeaway

The Jaccard Index measures the overlap between actual and predicted positive samples:

\[
\boxed{
J=\frac{TP}{TP+FP+FN}
}
\]

It is a useful metric when the quality of positive-class predictions matters more than the number of correctly classified negative samples.
