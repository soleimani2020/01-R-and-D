# Cosine Similarity Calculator

**Difficulty:** Easy  
**Topic:** Linear Algebra / Vectors

## Description

Write a Python function that calculates the **cosine similarity** between two vectors without using built-in functions like NumPy.  

- `v1` and `v2` are lists of integers or floats.  
- Cosine similarity measures the angle between two vectors and returns a value between **-1 and 1**.  
- If either vector is a zero vector, return **0.0**.



---

## Example

```python
from main import cosine_similarity

a = [1, 2, 3]
b = [4, 5, 6]
print(cosine_similarity(a, b))
# Output: 0.9746318461970762

# Zero vector example
c = [0, 0, 0]
print(cosine_similarity(a, c))
# Output: 0.0
