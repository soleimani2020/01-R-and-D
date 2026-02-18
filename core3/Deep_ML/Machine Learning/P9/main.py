def matrix_multiply(A: list[list[float]], B: list[list[float]]) -> list[list[float]] | int:
    """
    Multiply two matrices A and B.
    
    Args:
        A: Matrix of size m x n
        B: Matrix of size p x q
    
    Returns:
        The product matrix C = A * B if shapes align,
        otherwise returns -1
    """
    
    # Check if multiplication is possible
    if len(A[0]) != len(B):
        return -1  # number of columns in A != number of rows in B
    
    m, n = len(A), len(A[0])
    p, q = len(B), len(B[0])
    
    # Initialize result matrix with zeros
    C = [[0 for _ in range(q)] for _ in range(m)]
    
    # Multiply
    for i in range(m):
        for j in range(q):
            for k in range(n):  # or range(p), same as n
                C[i][j] += A[i][k] * B[k][j]
    
    return C
