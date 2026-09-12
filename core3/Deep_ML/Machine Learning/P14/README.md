# Calculate Covariance Matrix

**Difficulty:** Easy
**Topic:** Statistics

Write a Python function to calculate the **covariance matrix** for a set of features.

Each inner list represents one feature and its observations.

## Formula

For two variables \(X\) and \(Y\):

$$
\operatorname{Cov}(X,Y)=
\frac{1}{n-1}
\sum_{i=1}^{n}(x_i-\bar{x})(y_i-\bar{y})
$$


Output:

```python
[
    [1.0, 2.0],
    [2.0, 4.0]
]
```

## Test Cases

```python
assert covariance_matrix([[1, 2, 3], [2, 4, 6]]) == [
    [1.0, 2.0],
    [2.0, 4.0]
]
```
