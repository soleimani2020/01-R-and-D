# Calculate Jaccard Index for Binary Classification

## Overview

The **Jaccard Index** measures the similarity between two binary sets.

For binary classification, it compares the overlap between:

- `y_true` = true labels
- `y_pred` = predicted labels

The Jaccard Index ranges from:

$$
0 \le J \le 1
$$

where:

- `0` means no overlap
- `1` means perfect overlap

---

## Formula

The Jaccard Index is:

$$
J =
\frac{\text{Intersection}}{\text{Union}}
$$

For binary classification:

$$
J =
\frac{TP}{TP + FP + FN}
$$

where:

- $TP$ = True Positives
- $FP$ = False Positives
- $FN$ = False Negatives

True negatives are not included.

---
