# Matrix Multiplication

**Difficulty:** Medium  
**Topic:** Linear Algebra / Matrices

## Description

Write a Python function that **multiplies two matrices**.  

- `A` and `B` are 2D lists of numbers (integers or floats).  
- The function should return the **product matrix `C = A ⋅ B`** if the dimensions are compatible.  
- If the number of columns in `A` is **not equal** to the number of rows in `B`, return **-1**.  

Matrix multiplication rules:

- If `A` is an \(n × m\) matrix and `B` is an \(m × p\) matrix, the result `C` is an \(n × p\) matrix.  
- Each element of `C` is calculated as:

\[
C[i][j] = \sum_{k=0}^{m-1} A[i][k] \cdot B[k][j]
\]

---

## Example

```python
from main import matrix_multiply

A = [
    [1, 2],
    [3, 4]
]
B = [
    [5, 6],
    [7, 8]
]

C = matrix_multiply(A, B)
print(C)
# Output:
# [
#   [19, 22],
#   [43, 50]
# ]

# Incompatible shapes
D = [
    [1, 2, 3],
    [4, 5, 6]
]
E = [
    [7, 8],
    [9, 10]
]
result = matrix_multiply(D, E)
print(result)
# Output: -1

Video Tutorial: https://www.youtube.com/watch?v=LGiznhmpxp4
