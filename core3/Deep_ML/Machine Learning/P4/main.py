from typing import List, Union

def scalar_multiply(matrix: List[List[Union[int, float]]], scalar: Union[int, float]) -> List[List[Union[int, float]]]:
    """
    Multiply each element of a matrix by a scalar in-place.

    Parameters:
        matrix (list of list of int/float): Input matrix
        scalar (int/float): Scalar to multiply

    Returns:
        list of list of int/float: The same matrix after scalar multiplication
    """
    if not matrix or not matrix[0]:
        raise ValueError("Matrix cannot be empty")
    
    p, q = len(matrix), len(matrix[0])  # Rows / Columns
    
    for i in range(p):
        for j in range(q):
            matrix[i][j] = matrix[i][j] * scalar
    
    return matrix

# Example usage
if __name__ == "__main__":
    mat = [[1, 2, 3],
           [4, 5, 6],
           [7, 8, 9]]
    scalar = 2
    
    print("Original matrix:")
    for row in mat:
        print(row)
    
    scalar_multiply(mat, scalar)
    
    print("\nMatrix after multiplying by", scalar, ":")
    for row in mat:
        print(row)
