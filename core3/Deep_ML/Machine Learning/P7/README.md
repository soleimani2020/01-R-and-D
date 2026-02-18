# Eigenvalues of a 2x2 Matrix

**Difficulty:** Medium  
**Topic:** Linear Algebra / Matrices

## Description

Write a Python function that calculates the **eigenvalues** of a **2x2 matrix**.  

- `matrix`: a 2×2 list of lists representing the matrix  
- The function should return a **list of two eigenvalues**, sorted from **highest to lowest**.

### Formula for a 2x2 matrix

For a matrix:

\[
A = \begin{bmatrix} a & b \\ c & d \end{bmatrix}
\]

The eigenvalues are the roots of the characteristic equation:

\[
\lambda^2 - (a + d)\lambda + (ad - bc) = 0
\]

Solve for \(\lambda\) to get the two eigenvalues.

---

## Example

```python
from main import calculate_eigenvalues

matrix = [
    [4, 2],
    [1, 3]
]

eigenvalues = calculate_eigenvalues(matrix)
print(eigenvalues)
# Output: [5.0, 2.0]  # sorted from highest to lowest

