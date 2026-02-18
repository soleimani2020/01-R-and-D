import numpy as np

def calculate_dot_product(vec1, vec2):
	"""
	Calculate the dot product of two vectors.
	Args:
		vec1 (numpy.ndarray): 1D array representing the first vector.
		vec2 (numpy.ndarray): 1D array representing the second vector.
	Returns:
		The dot product of the two vectors.
	"""
	# Your code here
	
	
	if len(vec1) != len(vec2):
	    raise ValueError("Vectors must have the same length") 
	
	res = 0
	for i in range(len(vec1)):
	    res+=vec1[i]*vec2[i]
	    
	return res
	
	
	

# Example
v1 = np.array([1, 2, 3])
v2 = np.array([4, 5, 6])
print("Dot product:", calculate_dot_product(v1, v2))
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
