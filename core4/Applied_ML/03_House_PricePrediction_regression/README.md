# Project 03: House Price Prediction with Linear and Regularized Regression

> 🏡 **Machine Learning** · 📊 **Regression Analysis** · ⚙️ **Feature Engineering** · 🎯 **Model Evaluation** · 🛡️ **Regularization**

---

## 📌 Project Overview

House price prediction is a fundamental regression task in machine learning that involves estimating property values from characteristics such as living area, overall quality, construction year, and other structural and categorical attributes.

This project explores the Ames Housing dataset to predict residential sale prices using three regression algorithms: **Linear Regression, Ridge Regression, and Elastic Net**.

The workflow uses scikit-learn pipelines to integrate missing-value imputation, feature scaling, categorical encoding, cross-validation, and model selection into a reproducible and leakage-free process.

The project also investigates regularization, feature selection, prediction errors, residual patterns, and the limitations of applying a linear model to real-estate prices.

---

## 🎯 Main Goal

Develop and evaluate regression models capable of predicting house sale prices from property characteristics, comparing ordinary least squares with Ridge and Elastic Net regularization using cross-validation and a held-out test set.

---

## 🎯 Project Objectives

- 🏡 **Price Prediction:** Estimate residential sale prices using structural and categorical property features.
- 🔍 **Data Exploration:** Inspect the dataset, variable types, missing values, and relationships with sale price.
- 🧹 **Data Preprocessing:** Handle missing values and distinguish numerical measurements from categorical building codes.
- ⚙️ **Feature Engineering:** Build a preprocessing pipeline with imputation, standardization, and one-hot encoding.
- 📈 **Linear Regression:** Establish an ordinary least-squares regression benchmark.
- 🛡️ **Regularization:** Compare Ridge's L2 penalty with Elastic Net's L1/L2 regularization framework.
- 🔬 **Hyperparameter Tuning:** Use five-fold cross-validation to select regularization parameters.
- 📏 **Model Evaluation:** Compare CV RMSE, test MAE, test RMSE, and test R².
- 📊 **Error Analysis:** Examine actual versus predicted prices, residuals, and performance on expensive homes.
- 🎯 **Feature Selection:** Investigate the sparse coefficients produced by the selected Elastic Net model.

---

## 📊 Key Results

The selected Elastic Net model achieves a test R² of **0.898**, explaining approximately 90% of the variation in sale prices across the held-out test set.

| Model | CV RMSE | Test MAE | Test RMSE | Test R² |
|---|---:|---:|---:|---:|
| **Elastic Net** | **$34,075** | **$17,055** | **$26,623** | **0.898** |
| Ridge | $35,078 | $18,999 | $29,475 | 0.875 |
| Linear Regression | $39,208 | $18,450 | $28,102 | 0.887 |
| Mean baseline | $77,067 | $59,547 | $83,545 | 0.000 |

### 🏆 Best-Performing Model

**Selected model: Elastic Net**

- **Alpha:** `100`
- **L1 ratio:** `1.0`
- **CV RMSE:** `$34,075`
- **Test MAE:** `$17,055`
- **Test RMSE:** `$26,623`
- **Test R²:** `0.898`

The model was selected using cross-validated RMSE only. The test set was reserved for final evaluation and was not used to choose the model or its hyperparameters.

### 🔑 Key Findings

- **Prediction accuracy:** The median absolute error is approximately **$12,000**, while the test RMSE is approximately **$27,000**. RMSE is more sensitive to large prediction errors.
- **Feature selection:** The selected model retains **111 of 298 input columns** through its L1 penalty.
- **Main price drivers:** Above-ground living area (`GrLivArea`), overall quality (`OverallQual`), and year built (`YearBuilt`) are important predictive features.
- **Regularization:** Elastic Net outperforms Linear Regression in both cross-validation and test RMSE.
- **Ridge comparison:** Ridge performs better than Linear Regression in cross-validation but worse on the held-out test set. The difference should be interpreted cautiously because it may depend on the particular data split.

---

## 🧭 Methodology

| Step | Description | Purpose |
|---|---|---|
| 1. Data Loading | Load 1,460 homes and 79 original features | Establish the dataset for analysis |
| 2. Data Inspection | Examine variable types, missingness, and target distribution | Understand data quality and structure |
| 3. Train-Test Split | Allocate 70% to training (1,022 homes) and 30% to testing (438 homes) using a fixed random seed | Reserve unseen observations for final evaluation |
| 4. Feature Preparation | Convert `MSSubClass` from a numeric building code to a categorical feature | Avoid treating category codes as continuous measurements |
| 5. Missing-Value Handling | Drop columns with more than 80% missing values based on training data; impute remaining values | Reduce missing-data problems without using test information |
| 6. Preprocessing Pipeline | Apply median imputation, missing indicators, and standardization to numeric features; use `"None"` imputation and one-hot encoding for categorical features | Prepare features consistently |
| 7. Baseline Model | Predict the mean sale price | Establish a reference point for regression performance |
| 8. Model Training | Fit Linear Regression, Ridge, and Elastic Net | Compare unregularized and regularized regression |
| 9. Hyperparameter Tuning | Apply five-fold `GridSearchCV` for Ridge alpha and Elastic Net alpha and L1 ratio | Select hyperparameters using training data only |
| 10. Model Selection | Choose the model with the lowest mean cross-validated RMSE | Prevent test-set information from influencing selection |
| 11. Final Evaluation | Calculate test MAE, RMSE, and R² | Measure performance on unseen homes |
| 12. Diagnostics | Compare training and test errors; plot predictions and residuals; inspect coefficients | Investigate generalization, error patterns, and feature importance |

All learned preprocessing steps are contained within the scikit-learn pipeline, allowing preprocessing to be refitted within each cross-validation fold.

---

## 🧠 Models and Regularization

### 1. Linear Regression

Linear Regression estimates coefficients by minimizing the sum of squared prediction errors.

It provides a baseline for understanding how well a linear combination of property characteristics can explain sale prices without an explicit coefficient penalty.

### 2. Ridge Regression

Ridge Regression adds an L2 penalty to the least-squares objective.

This penalty shrinks coefficients toward zero, helping control coefficient magnitude and reduce sensitivity to correlated predictors.

### 3. Elastic Net

Elastic Net combines L1 and L2 regularization. Its L1 component can drive coefficients to exactly zero, while its L2 component encourages coefficient shrinkage.

The selected model uses `alpha = 100` and `l1_ratio = 1.0`. At this L1 ratio, Elastic Net reduces mathematically to **Lasso-style L1 regularization** rather than a mixture of L1 and L2 penalties.

The selected model retains **111 of 298 input columns**, demonstrating its ability to produce a sparse representation of the encoded features.

| Model | Regularization | Main purpose |
|---|---|---|
| Linear Regression | None | Establish an ordinary least-squares baseline |
| Ridge | L2 | Shrink coefficients and control model complexity |
| Elastic Net (`l1_ratio = 1.0`) | L1 | Perform coefficient shrinkage and feature selection |

---

## 🗂️ Dataset

The project uses the **Ames Housing dataset**, documented by De Cock (2011) and popularized through Kaggle's *House Prices: Advanced Regression Techniques* competition.

The dataset contains **1,460 residential properties and 79 original features**, with `SalePrice` as the prediction target.

Features cover several aspects of each property, including:

- **Building characteristics:** Building type, construction year, and structural attributes.
- **Size and layout:** Above-ground living area, room configuration, and basement-related measurements.
- **Quality and condition:** Overall material and finish quality, general condition, and related attributes.
- **Location and classification:** Neighborhood and property classification variables.
- **Amenities:** Garage, fireplace, porch, and other property features.

### Data Preparation

- Convert `MSSubClass` to a categorical feature.
- Determine columns exceeding the 80% missingness threshold using training data only.
- Impute numerical missing values using the median and add missing-value indicators.
- Fill categorical missing values with an explicit `"None"` category where missingness represents the absence of a feature.
- Apply standardization to numerical features.
- Apply one-hot encoding to categorical features.
- Keep the held-out test set separate from model selection and cross-validation.

The original dataset is split into **1,022 training homes (70%)** and **438 test homes (30%)**.

---

## 🛡️ Leakage Prevention and Generalization

Data leakage occurs when information from validation or test observations influences model training, feature preparation, or model selection.

This project reduces leakage risk through the following practices:

- **Split first:** Separate the training and test observations before fitting preprocessing steps.
- **Pipeline-based preprocessing:** Place imputation, scaling, and encoding inside the scikit-learn pipeline.
- **Fold-specific transformations:** Refit preprocessing within each cross-validation fold.
- **Training-only selection:** Use cross-validated RMSE to select the model and hyperparameters.
- **Held-out final evaluation:** Use the test set only after model selection.

These practices help provide a more reliable estimate of performance on unseen observations.

---

## 📏 Model Evaluation

The project evaluates models using complementary metrics that measure different aspects of prediction quality.

| Metric | Description |
|---|---|
| **CV RMSE** | Root mean squared error averaged across cross-validation folds; used for model selection |
| **Test MAE** | Average absolute prediction error on the held-out test set |
| **Test RMSE** | Root mean squared prediction error, penalizing large errors more heavily |
| **Test R²** | Proportion of test-target variation explained relative to a mean-based reference |

### Results Interpretation

**Elastic Net** achieves the lowest cross-validation RMSE ($34,075) and test RMSE ($26,623) among the three fitted regression models.

Its test R² of `0.898` indicates strong predictive performance on this particular test split. However, R² does not mean that 89.8% of individual predictions are accurate, nor does it establish causal relationships between property characteristics and prices.

The mean baseline performs substantially worse, with a test MAE of $59,547 and a test RMSE of $83,545.

---

## 🔬 Model Interpretation and Error Analysis

Understanding prediction errors is important because aggregate metrics can hide systematic weaknesses.

### Important Predictive Features

The analysis identifies the following as major price-related predictors:

- **`GrLivArea`:** Above-ground living area.
- **`OverallQual`:** Overall material and finish quality.
- **`YearBuilt`:** Year of construction.

These features provide useful predictive information, but their coefficients should not be interpreted as causal effects.

### Residual Diagnostics

The diagnostic analysis identifies two important patterns:

- **Under-prediction of expensive homes:** Homes in the top 10% of predicted prices are under-predicted by approximately **$8,000 on average**, according to the reported analysis.
- **Increasing residual spread:** Errors fan out as prices rise, suggesting that error variance is not constant across the price range.

The median absolute error of approximately $12,000 is lower than the RMSE of approximately $27,000, indicating that larger errors have a noticeable effect on the squared-error metric.

### Rare Categories

Some one-hot encoded categories occur very infrequently. For example, `RoofMatl_ClyTile` occurs in only one training home.

Coefficients associated with such rare categories can be unstable and should not be treated as reliable general rules about house prices.

---

## ⚠️ Limitations

- **Single test split:** The 438-home test set provides only one assessment of generalization. Small performance differences may vary across train-test splits.
- **Limited geographic scope:** The data represent Ames, Iowa, over the 2006–2010 period. The results should not be assumed to generalize to other cities or current property markets without further validation.
- **Linear modeling assumptions:** Modeling raw sale price with linear regression can struggle with nonlinear relationships and increasing residual variance.
- **Expensive-home bias:** The model under-predicts expensive homes in the reported diagnostics.
- **Rare-category instability:** Coefficients for categories represented by very few observations may be unreliable.
- **Feature selection trade-offs:** Sparse coefficients simplify the model, but correlated features and the selected regularization strength can affect which variables remain nonzero.

---

## 🚀 Next Steps

Several improvements could strengthen the model and provide more robust comparisons.

- **Log-transform the target:** Model `log(SalePrice)` to investigate whether the transformation reduces the fan-shaped residual pattern and improves predictions for expensive properties.
- **Repeated cross-validation:** Evaluate model performance across repeated splits to obtain more stable comparisons.
- **Pool rare categories:** Group infrequent categorical levels to reduce the instability of their coefficients.
- **Explore nonlinear models:** Compare the linear approaches with tree-based regression methods to capture nonlinear relationships and interactions.
- **Evaluate residual subgroups:** Compare errors across price ranges, neighborhoods, and property-quality groups.
- **Assess practical error:** Examine MAE and prediction errors by price band to understand where the model is most useful.

Any changes should be assessed using training-only cross-validation, with the test set reserved for final evaluation.

---

## 🛠️ Technologies & Tools

| Category | Technologies |
|---|---|
| Programming Language | Python |
| Data Manipulation | Pandas, NumPy |
| Machine Learning | Scikit-learn |
| Regression Models | Linear Regression, Ridge, Elastic Net |
| Preprocessing | `Pipeline`, `ColumnTransformer`, `SimpleImputer`, `OneHotEncoder`, `StandardScaler` |
| Hyperparameter Tuning | `GridSearchCV`, five-fold cross-validation |
| Model Evaluation | MAE, RMSE, R² |
| Visualization | Matplotlib |
| Development Environment | Jupyter Notebook |

---

## 📁 Project Structure

```text
03_House_Price_Prediction/
│
├── README.md
├── House_Price.ipynb
├── house_price.csv
├── house-price-project.png
│
├── figures/                 # Prediction, residual, and diagnostic plots
├── models/                  # Saved model artifacts, if included
└── results/                 # Metrics and analysis outputs, if included
```

The structure can be adapted to match the actual repository. Large datasets and generated model artifacts do not necessarily need to be committed to version control.

---

## 🚀 How to Run

Clone the repository and install the required dependencies:

```bash
git clone https://github.com/<your-username>/<your-repo>.git
cd <your-repo>

pip install numpy pandas matplotlib scikit-learn jupyter

jupyter notebook House_Price.ipynb
```

The notebook reads `house_price.csv` when it is available in the same directory; otherwise, it downloads the dataset from the configured GitHub source.

Use **scikit-learn 1.2 or newer**, because the notebook uses the `OneHotEncoder(sparse_output=...)` parameter.

---

## 🔬 Scientific Motivation

House price prediction illustrates how supervised learning can transform structured property data into quantitative estimates.

Linear regression provides an interpretable baseline, while regularization helps manage the complexity introduced by correlated numerical variables and high-dimensional one-hot encoded categorical features.

By combining leakage-aware preprocessing, cross-validation, held-out testing, and residual diagnostics, this project demonstrates an end-to-end workflow for building and evaluating regression models on real-world tabular data.

---

## 📌 Project Focus

**Ames Housing → Exploratory Analysis → Leakage-Free Preprocessing → Linear Regression → Ridge & Elastic Net → Cross-Validation → Model Selection → Test Evaluation → Residual Analysis → Feature Interpretation**

---

## 👤 Author

**<Your name>** · [LinkedIn](https://www.linkedin.com/in/your-profile) · [Portfolio](https://your-site.example)
