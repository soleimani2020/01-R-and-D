from collections import Counter

def confusion_matrix(data):
    counts = Counter((y_true, y_pred) for y_true, y_pred in data)
    
    tp = counts[(1, 1)]
    tn = counts[(0, 0)]
    fp = counts[(0, 1)]
    fn = counts[(1, 0)]
    
    return [[tp, fn],
            [fp, tn]]


def confusion_matrix(data):
    tp = tn = fp = fn = 0
    
    for y_true, y_pred in data:
        if y_true == 1 and y_pred == 1:
            tp += 1
        elif y_true == 1 and y_pred == 0:
            fn += 1
        elif y_true == 0 and y_pred == 1:
            fp += 1
        elif y_true == 0 and y_pred == 0:
            tn += 1
    
    return [[tp, fn],
            [fp, tn]]
