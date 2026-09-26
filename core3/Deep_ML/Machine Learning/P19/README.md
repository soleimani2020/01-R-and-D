# Generate a Confusion Matrix for Binary Classification

**Difficulty:** Easy  
**Topic:** Machine Learning

## Task

Implement the function:

```python
confusion_matrix(data)
```

The function should generate a **2 × 2 confusion matrix** for a binary classification problem.

Each observation contains:

```python
[y_true, y_pred]
```

where:

- `y_true` is the actual label.
- `y_pred` is the predicted label.
- Labels are either `0` or `1`.

## Confusion Matrix

For binary classification, the confusion matrix has the following structure:

```text
                Predicted
                0       1
Actual  0      TN      FP
        1      FN      TP
```

where:

- **TN — True Negative:** actual `0`, predicted `0`
- **FP — False Positive:** actual `0`, predicted `1`
- **FN — False Negative:** actual `1`, predicted `0`
- **TP — True Positive:** actual `1`, predicted `1`

The function should therefore return:

```python
[
    [TN, FP],
    [FN, TP]
]
```

## Input

A list of lists, where each inner list contains:

```python
[y_true, y_pred]
```

Example:

```python
data = [
    [1, 1],
    [0, 0],
    [1, 0],
    [0, 1],
    [1, 1]
]
```

## Output

A `2 × 2` confusion matrix represented as a list of lists.

For the example above:

```python
[
    [1, 1],
    [1, 2]
]
```

Explanation:

- `TN = 1`
- `FP = 1`
- `FN = 1`
- `TP = 2`

