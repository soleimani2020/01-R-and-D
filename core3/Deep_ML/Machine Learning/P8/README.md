# 2x2 Matrix Inverse Calculator

**Difficulty:** Easy  
**Topic:** Linear Algebra / Matrices

## Description

Write a Python function that calculates the **inverse of a 2x2 matrix**.  

- `matrix`: a 2×2 list of lists representing the matrix  
- The inverse of a matrix \(A\) is another matrix \(A_{\text{inv}}\) such that:

\[
A \cdot A_{\text{inv}} = I
\]

where \(I\) is the 2x2 identity matrix.

- For a 2x2 matrix:

\[
A = \begin{bmatrix} a & b \\ c & d \end{bmatrix}
\]

the inverse exists **only if the determinant** \(ad - bc \neq 0\).  
- If the determinant is zero, return **None**.

---

## Example

```python
from main import inverse_2x2_matrix

matrix = [
    [4, 7],
    [2, 6]
]

inverse = inverse_2x2_matrix(matrix)
print(inverse)
# Output:
# [
#   [0.6, -0.7],
#   [-0.2, 0.4]
# ]

# Non-invertible matrix
matrix_singular = [
    [1, 2],
    [2, 4]
]
inverse = inverse_2x2_matrix(matrix_singular)
print(inverse)
# Output: None


Video Tutorial: https://www.youtube.com/watch?v=-fhFySMHPZk
