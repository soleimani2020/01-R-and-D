def calculate_matrix_mean(matrix: list[list[float]], mode: str) -> list[float]:
    """
    Calculate the mean of rows or columns of a matrix.

    Args:
        matrix (list[list[float]]): 2D matrix
        mode (str): 'row' for row-wise mean, 'col' for column-wise mean

    Returns:
        list[float]: List of means for each row or column
    """
    
    rows, cols = len(matrix), len(matrix[0])
    means = []
    
    if mode == 'row':
        for row in matrix:
            mean = sum(row) / len(row)
            means.append(mean)
    elif mode == 'col':
        for i in range(cols):
            column = [row[i] for row in matrix]
            mean = sum(column) / len(column)
            means.append(mean)
    else:
        raise ValueError("Mode must be 'row' or 'col'")
    
    return means
