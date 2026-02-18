def cosine_similarity(v1, v2):
    """
    Calculate the cosine similarity of two vectors without using built-in functions.

    Args:
        v1 (list of int/float): First vector
        v2 (list of int/float): Second vector

    Returns:
        float: Cosine similarity between -1 and 1
    """
    
    if len(v1) != len(v2):
        raise ValueError("Vectors must be of the same length")
    
    # Compute dot product : numerator = np.sum(v1 * v2)
    dot = 0
    for i in range(len(v1)):
        dot += v1[i] * v2[i]
    
    # Compute norm of v1 : np.linalg.norm(v1)
    norm1 = 0
    for i in range(len(v1)):
        norm1 += v1[i] * v1[i]
    norm1 = norm1 ** 0.5  # square root
    
    # Compute norm of v2
    norm2 = 0
    for i in range(len(v2)):
        norm2 += v2[i] * v2[i]
    norm2 = norm2 ** 0.5  # square root
    
    # Avoid division by zero
    if norm1 == 0 or norm2 == 0:
        return 0.0
    
    # Cosine similarity
    return dot / (norm1 * norm2)
