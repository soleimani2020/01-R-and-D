import numpy as np

def accuracy_score(y_true, y_pred):
	total_predictions = len(y_true)


	y_true = np.array(y_true)
	y_pred = np.array(y_pred)

	if y_true.shape != y_pred.shape:
		raise ValueError("y_true and y_pred must have the same shape")

	#correct_predictions = 0
	#for i in range(len(y_true)):
	#	if y_true[i] == y_pred[i]:
	#		correct_predictions+=1

	correct_predictions =np.sum(y_true==y_pred)

	return correct_predictions/total_predictions



	
