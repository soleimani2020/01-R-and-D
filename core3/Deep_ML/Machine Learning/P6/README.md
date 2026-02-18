# Matrix Mean Calculator

**Difficulty:** Easy  
**Topic:** Linear Algebra / Matrices

## Description

Write a Python function that calculates the **mean of each row or column** of a 2D matrix.  

- `matrix`: a 2D list of numbers (floats or integers)  
- `mode`: a string specifying the calculation mode:  
  - `'row'` → calculate the mean of each row  
  - `'col'` → calculate the mean of each column  

The function should return a **list of means** corresponding to each row or column.

---

## Example

```python
from main import calculate_matrix_mean

matrix = [
    [1.0, 2.0, 3.0],
    [4.0, 5.0, 6.0],
    [7.0, 8.0, 9.0]
]

row_means = calculate_matrix_mean(matrix, 'row')
print(row_means)
# Output: [2.0, 5.0, 8.0]

col_means = calculate_matrix_mean(matrix, 'col')
print(col_means)
# Output: [4.0, 5.0, 6.0]


Video Tutorial: https://www.youtube.com/watch?v=5KCS4qVczJ0

