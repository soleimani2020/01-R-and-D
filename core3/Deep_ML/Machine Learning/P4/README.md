# Scalar Multiply

**Difficulty:** Easy  
**Topic:** Linear Algebra / Matrices

## Description

Write a Python function to multiply each element of a matrix by a scalar.  

The function should take:

- `matrix`: a 2D list of integers or floats  
- `scalar`: an integer or float  

and return the same matrix after multiplying each element by the scalar.  

## Example

```python
from main import scalar_multiply

matrix = [
    [1, 2, 3],
    [4, 5, 6],
    [7, 8, 9]
]
scalar = 2

scalar_multiply(matrix, scalar)

print(matrix)
# Output:
# [
#   [2, 4, 6],
#   [8, 10, 12],
#   [14, 16, 18]
# ]

