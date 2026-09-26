import numpy as np

def f_score(y_true, y_pred, beta):
    TP = np.sum((y_true == 1) & (y_pred == 1))
    FP = np.sum((y_true == 0) & (y_pred == 1))
    FN = np.sum((y_true == 1) & (y_pred == 0))

    if TP + FP == 0:
        precision = 0.0
    else:
        precision = TP / (TP + FP)

    if TP + FN == 0:
        recall = 0.0
    else:
        recall = TP / (TP + FN)

    denominator = (beta ** 2 * precision) + recall

    if denominator == 0:
        return 0.0

    score = (1 + beta ** 2) * (precision * recall) / denominator

    return round(score, 3)
