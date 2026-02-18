# Matrix-Vector Dot Product

**Difficulty:** Easy  
**Topic:** Linear Algebra

## Description

Write a Python function that computes the **dot product** of a matrix and a vector.  

- `matrix`: a 2D list of numbers (integers or floats)  
- `vector`: a 1D list of numbers (integers or floats)  

The function should return:

- A **list** representing the resulting vector if the operation is valid  
- `-1` if the matrix and vector dimensions are incompatible  

**Rules:**

- A matrix (n × m) can be dotted with a vector only if the vector has **length m**.  
- The resulting vector will have **length n**.

---

## Example

```python
from main import matrix_vector_dot

matrix = [
    [1, 2, 3],
    [4, 5, 6]
]
vector = [7, 8, 9]

result = matrix_vector_dot(matrix, vector)
print(result)
# Output: [50, 122]

# Incompatible dimensions
bad_vector = [1, 2]
result = matrix_vector_dot(matrix, bad_vector)
print(result)
# Output: -1

