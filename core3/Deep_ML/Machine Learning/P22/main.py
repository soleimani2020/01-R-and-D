import numpy as np

def recall(y_true, y_pred):
    """
    Calculate the recall metric for binary classification.
    """
    TP = 0 
    FN = 0
    
    for true, pred in zip(y_true, y_pred):
        if true == 1 and pred == 1:
            TP += 1 
        elif true == 1 and pred == 0:
            FN += 1 
            
    if TP + FN == 0:
        return 0.0
            
    return TP / (TP + FN)
