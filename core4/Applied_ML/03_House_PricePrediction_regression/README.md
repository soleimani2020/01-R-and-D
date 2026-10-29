# House Price Prediction with Linear and Regularized Regression

Predicting home sale prices on the Ames Housing dataset and comparing **Linear Regression**, **Ridge** and **Elastic Net** with a leakage-free scikit-learn pipeline.

![Project overview](house-price-project.png)

## Question

Can property characteristics be used to estimate a home's sale price, and how do ordinary linear regression, Ridge regression and Elastic Net compare on test data held out from model development?

## Key results

Property characteristics explain about **90% of the variation** in sale prices on held-out homes.

| Model | CV RMSE | Test MAE | Test RMSE | Test R² |
|---|---:|---:|---:|---:|
| **Elastic Net** | **$34,075** | **$17,055** | **$26,623** | **0.898** |
| Ridge | $35,078 | $18,999 | $29,475 | 0.875 |
| Linear Regression | $39,208 | $18,450 | $28,102 | 0.887 |
| Mean baseline | $77,067 | $59,547 | $83,545 | 0.000 |

- **Selected model:** Elastic Net (`alpha = 100`, `l1_ratio = 1.0`), chosen by cross-validated RMSE only. The test set was used once, for the final evaluation.
- **Typical error:** median absolute error of about **$12k** on the test set (RMSE is higher, about $27k, because a few large misses count extra).
- **Feature selection:** the L1 penalty keeps 111 of 298 input columns.
- **Main drivers:** above-ground living area (`GrLivArea`), overall quality (`OverallQual`) and year built (`YearBuilt`).
- **Regularization:** Elastic Net beats Linear Regression on both CV and test RMSE. Ridge beats Linear Regression in CV but not on the test set, so its advantage is small and within split-to-split noise.

## Approach

1. **Load and inspect** the data (1,460 homes, 79 features). `MSSubClass` is converted to a categorical column because it is a building-type code, not a quantity.
2. **Split first:** 70% train (1,022 homes) / 30% test (438 homes), fixed random seed.
3. **Missing values:** columns with more than 80% missing values (decided on the training data only) are dropped. For most other categorical columns, a missing value means "the house has no such feature", so gaps are filled with an explicit `"None"` category.
4. **One `Pipeline` for everything:** median imputation with missing-indicators and standardization for numeric columns; `"None"` fill and one-hot encoding for categorical columns. Because it lives inside the pipeline, preprocessing is re-fitted on every CV fold and no validation or test row influences it.
5. **Baseline and tuning:** a mean-predictor baseline, plain Linear Regression, and `GridSearchCV` (5-fold) for Ridge (`alpha`) and Elastic Net (`alpha`, `l1_ratio`).
6. **Selection:** lowest mean CV RMSE. The test set plays no part in this choice.
7. **Diagnostics:** train vs. test errors, actual vs. predicted and residual plots, and a coefficient analysis.

## Limitations

- The test set has 438 homes and is a single split; small differences between models are within split-to-split noise. Repeated CV would give tighter estimates.
- The model is linear on raw price: it under-predicts expensive homes (the homes ranked in the top 10% of predicted price are under-predicted by about $8k on average) and the residuals fan out as price rises.
- One-hot coefficients for very rare categories (for example `RoofMatl_ClyTile`, 1 training home) are unreliable and should not be read as general rules.
- The data cover one city and period (Ames, Iowa, 2006-2010); results should not be applied to other markets without re-validation.

## Next steps

- Model `log(SalePrice)` to reduce the fan-shaped errors and the under-prediction of expensive homes.
- Pool rare categories and use repeated cross-validation for tighter model comparisons.

## Repository structure

```
.
├── House_Price.ipynb       # full analysis (rename to match your file)
├── house_price.csv         # Ames Housing training data (1,460 rows)
├── house-price-project.png # overview image
└── README.md
```

## How to run

```bash
git clone https://github.com/<your-username>/<your-repo>.git
cd <your-repo>
pip install numpy pandas matplotlib scikit-learn jupyter
jupyter notebook House_Price.ipynb
```

The notebook reads `house_price.csv` if it is in the same folder; otherwise it downloads the same data from GitHub. Use scikit-learn 1.2 or newer (the notebook uses `OneHotEncoder(sparse_output=...)`).

## Tech stack

Python, pandas, NumPy, scikit-learn (`Pipeline`, `ColumnTransformer`, `GridSearchCV`), Matplotlib.

## Data

Ames Housing dataset (De Cock, 2011), as used in the Kaggle competition *House Prices: Advanced Regression Techniques*.

## Author

**<Your name>** · [LinkedIn](https://www.linkedin.com/in/your-profile) · [Portfolio](https://your-site.example)
