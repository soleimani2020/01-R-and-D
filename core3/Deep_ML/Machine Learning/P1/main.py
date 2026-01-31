def matrix_dot_vector(a: list[list[int|float]], b: list[int|float]) -> list[int|float]:
    # Return a list where each element is the dot product of a row of 'a' with 'b'.
    # If the number of columns in 'a' does not match the length of 'b', return -1.
    
    if len(a) == 0:
        return []
    
    n = len(a)      # number of rows
    m = len(a[0])   # number of columns
    
    if len(b) != m:
        return -1
    
    result = []
    for row in a:
        dot = 0
        for matrix_element, vector_element in zip(row, b):
            dot += matrix_element * vector_element
        result.append(dot)
    
    return result



matrix = [
    [1, 2, 3],
    [4, 5, 6]
]
vector = [7, 8, 9]

print(matrix_dot_vector(matrix, vector))
