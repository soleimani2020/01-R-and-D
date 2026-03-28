def transpose_matrix(a: list[list[int|float]]) -> list[list[int|float]]:
    """
    Transpose a 2D matrix by swapping rows and columns.
   
    Args:
        a: A 2D matrix of shape (m, n)
   
    Returns:
        The transposed matrix of shape (n, m)
    """
    # Your code here
    m , n = len(a) , len(a[0])  # number of rows and columns in the matrix [m*n]
    transposed_a = [ [0]*m for _ in range(n) ] # [n*m]
   
    for i in range(n):
        for j in range(m):
            transposed_a[i][j] = a[j][i]
   
    return transposed_a
   
   
    import numpy as np
    a = np.array(a)
    return a.T
