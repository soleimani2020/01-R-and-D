# Calculate Performance Metrics for a Classification Model

A Python implementation of several important **binary classification performance metrics**.

## Problem

Write a function:

```python
performance_metrics(actual, predicted)
```

that calculates:

- Confusion Matrix
- Accuracy
- F1 Score
- Specificity
- Negative Predictive Value

The input lists contain binary labels:

- `1` = positive class
- `0` = negative class

Both input lists must have the same length.

## Confusion Matrix

For binary classification:

```text
                 Predicted
               0           1
Actual  0      TN          FP
        1      FN          TP
```

So the confusion matrix is:

```python
[
    [TN, FP],
    [FN, TP]
]
```

where:

- **TP** = True Positive
- **TN** = True Negative
- **FP** = False Positive
- **FN** = False Negative

## Accuracy

$$
\text{Accuracy} =
\frac{TP + TN}
{TP + TN + FP + FN}
$$

Accuracy measures the proportion of all predictions that are correct.

## Precision

$$
\text{Precision} =
\frac{TP}{TP + FP}
$$

Precision is needed to calculate the F1 Score.

## Recall

$$
\text{Recall} =
\frac{TP}{TP + FN}
$$

Recall measures how many actual positive samples were correctly identified.

## F1 Score

$$
F_1 =
2 \times
\frac{\text{Precision} \times \text{Recall}}
{\text{Precision} + \text{Recall}}
$$

An equivalent formula using TP, FP, and FN is:

$$
F_1 =
\frac{2TP}
{2TP + FP + FN}
$$

## Specificity

$$
\text{Specificity} =
\frac{TN}{TN + FP}
$$

Specificity measures how effectively the model identifies actual negative samples.

## Negative Predictive Value

$$
NPV =
\frac{TN}{TN + FN}
$$

Negative Predictive Value measures how often a negative prediction is actually correct.


