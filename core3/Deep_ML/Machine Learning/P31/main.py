import numpy as np

def jaccard_index(y_true, y_pred):
    y_true = np.array(y_true)
    y_pred = np.array(y_pred)

    intersection = np.sum((y_true == 1) & (y_pred == 1))
    union = np.sum((y_true == 1) | (y_pred == 1))

    if union == 0:
        return 1.0

    result = intersection / union

    return round(float(result), 3)
