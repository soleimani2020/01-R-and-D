# Advanced Loan Classification and Risk Analysis

## Overview
This project develops and optimizes machine learning classification models to predict loan repayment outcomes using historical lending data from 2015 and 2016.

The objective is to identify loans with a higher likelihood of repayment, with particular attention to loans offering higher interest rates. The project emphasizes appropriate evaluation metrics, model optimization, and the prevention of data leakage.

## Objectives
- Explore and preprocess historical lending data.
- Identify variables available at the time of loan application.
- Prevent data leakage from information recorded after loan issuance.
- Build and compare classification models.
- Optimize model hyperparameters.
- Evaluate performance using suitable classification metrics.
- Investigate the trade-off between potential returns and repayment risk.

## Methods
- Exploratory data analysis
- Data cleaning and feature engineering
- Supervised machine learning
- Classification and model comparison
- Cross-validation and hyperparameter tuning
- Precision, recall, F1-score, and ROC-AUC evaluation
- Confusion matrix analysis
- Data leakage prevention
- Risk-oriented model evaluation

## Technologies
- Python
- Pandas
- NumPy
- Matplotlib
- Seaborn
- Scikit-learn

## Dataset
The project uses `lending_file_2015_2016.csv`, which contains information about loans requested and accepted during 2015 and 2016.

The dataset includes variables recorded at application time and others added later. Careful feature selection is therefore necessary to avoid using information that would not have been available when making a lending decision.

## Project Structure
```text
10_Advanced_Loan_Classification/
├── README.md
└── project.ipynb
```

## Status
Work in progress.
