import numpy as np

def predict_logistic(X: np.ndarray, weights: np.ndarray, bias: float) -> np.ndarray:
	"""
	Implements binary classification prediction using Logistic Regression.

	Args:
		X: Input feature matrix (shape: N x D)
		weights: Model weights (shape: D)
		bias: Model bias

	Returns:
		Binary predictions (0 or 1)
	"""
	# Linear score
	z = X @ weights + bias 
	
	# Sigmoid probability
	sigma = 1 / (1 + np.exp(-z))
	
	#np.where(condition, value_if_true, value_if_false)
	ans = np.where(sigma >= 0.5, 1, 0)

	return ans 
	 
	
