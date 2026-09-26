# Implement K-Fold Cross-Validation

A simple NumPy implementation of **K-Fold Cross-Validation** index splitting.

## Problem

Write a function that divides dataset indices into `k` folds and returns a list of train-test index pairs.

### Inputs

- `n_samples`: total number of samples
- `k`: number of folds, default `5`
- `shuffle`: whether to shuffle indices before splitting, default `True`

### Output

A list of `k` tuples:

```python
(train_indices, test_indices)
```

where both are Python lists of integers.

## Requirements

- Split the indices into `k` roughly equal folds.
- If `n_samples` is not divisible by `k`, distribute the extra samples to the first folds.
- Use one fold as the test set and all remaining folds as the training set.
- If `shuffle=True`, shuffle the indices with `np.random.shuffle()`.
- If `shuffle=False`, keep the original order.

## Example Fold Sizes

For:

```text
n_samples = 10
k = 3
```

the fold sizes should be:

```text
[4, 3, 3]
```

because:

```text
10 // 3 = 3
10 % 3  = 1
```

So the first fold receives one extra sample.

