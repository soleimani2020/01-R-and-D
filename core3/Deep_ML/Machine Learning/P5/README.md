# Cosine Similarity Calculator

**Difficulty:** Easy  
**Topic:** Linear Algebra / Vectors

## Description

Write a Python function that calculates the **cosine similarity** between two vectors without using built-in functions like NumPy.  

- `v1` and `v2` are lists of integers or floats.  
- Cosine similarity measures the angle between two vectors and returns a value between **-1 and 1**.  
- If either vector is a zero vector, return **0.0**.

### Cosine Similarity Formula

\[
\text{cosine\_similarity} = \frac{v1 \cdot v2}{||v1|| \cdot ||v2||}
\]

where \(v1 \cdot v2\) is the dot product, and \(||v1||\), \(||v2||\) are vector norms.

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
