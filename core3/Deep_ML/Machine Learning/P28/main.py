def performance_metrics(actual, predicted):
    TP = 0
    TN = 0
    FP = 0
    FN = 0

    for true, pred in zip(actual, predicted):
        if true == 1 and pred == 1:
            TP += 1
        elif true == 0 and pred == 0:
            TN += 1
        elif true == 0 and pred == 1:
            FP += 1
        elif true == 1 and pred == 0:
            FN += 1

    confusion_matrix = [
        [TN, FP],
        [FN, TP]
    ]

    total = TP + TN + FP + FN

    accuracy = (TP + TN) / total if total != 0 else 0.0

    f1_score = (
        (2 * TP) / (2 * TP + FP + FN)
        if (2 * TP + FP + FN) != 0
        else 0.0
    )

    specificity = (
        TN / (TN + FP)
        if (TN + FP) != 0
        else 0.0
    )

    negative_predictive_value = (
        TN / (TN + FN)
        if (TN + FN) != 0
        else 0.0
    )

    return (
        confusion_matrix,
        accuracy,
        f1_score,
        specificity,
        negative_predictive_value
    )
