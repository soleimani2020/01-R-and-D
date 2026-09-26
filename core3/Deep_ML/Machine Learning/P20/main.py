import numpy as np

def accuracy_score(y_true, y_pred):
	# Your code here
	correct = np.sum(y_true == y_pred)
	Total = len(y_true)
	
	accuracy = correct / Total
	
	return accuracy
	    
	    
	    
def accuracy_score(y_true, y_pred):
    correct = 0

    for true, pred in zip(y_true, y_pred):
        if true == pred:
            correct += 1

    return correct / len(y_true)	    
	    
	    
	    
y_true = np.array([1, 0, 1, 1, 0, 1])
y_pred = np.array([1, 0, 0, 1, 0, 1])
output = accuracy_score(y_true, y_pred)
print(output)
