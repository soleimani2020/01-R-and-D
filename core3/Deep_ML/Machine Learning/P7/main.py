import math

def calculate_eigenvalues(matrix: list[list[float|int]]) -> list[float]:
    a, b = matrix[0]  # First row
    c, d = matrix[1]  # Second row
    
    trace = a + d
    det = a*d - b*c
    
    discriminant = trace**2 - 4*det
    
    if discriminant < 0:
        raise ValueError("Matrix has complex eigenvalues, only real eigenvalues are supported")
    
    sqrt_discriminant = math.sqrt(discriminant)
    
    lambda1 = (trace + sqrt_discriminant) / 2
    lambda2 = (trace - sqrt_discriminant) / 2
    
    return sorted([lambda1, lambda2], reverse=True)
