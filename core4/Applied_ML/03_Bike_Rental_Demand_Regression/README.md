# Project 03: Bike Rental Demand Prediction with Regression

> 🚲 **Supervised Learning** · 📊 **Regression Analysis** · ⚙️ **Regularization** · 📈 **Predictive Modeling**

---

## 📌 Project Overview

This project applies **supervised machine learning regression techniques** to predict hourly bike rental demand using calendar, time, and weather-related features.

Using the **UCI Bike Sharing hourly dataset**, the project follows an end-to-end machine learning workflow, including data inspection, exploratory data analysis, chronological train/test splitting, preprocessing, regression modeling, hyperparameter tuning, and performance evaluation.

The models learn relationships between factors such as hour of the day, temperature, and other recorded conditions and the total number of bikes rented in each hour.

### 🎯 Main Goal

Build and evaluate regression models that predict hourly bike rental demand, with attention to predictive accuracy, regularization, chronological evaluation, and the interpretation of prediction errors.

**Target variable:** `cnt` — the total number of bike rentals in an hour.

---

## 🎯 Project Objectives

- 🚲 **Dataset Exploration** — Inspect the bike-sharing dataset, its features, data types, and target variable.
- 🔎 **Exploratory Data Analysis** — Explore rental patterns by hour and relationships between weather-related features and demand.
- 🧹 **Data Preprocessing** — Check missing values and prepare numerical and categorical features.
- ⚙️ **Feature Engineering** — One-hot encode categorical features and standardize numerical features.
- 🕒 **Chronological Evaluation** — Preserve time order and reserve the final 20% of observations for testing.
- 📊 **Model Development** — Train Linear Regression, Ridge Regression, and Elastic Net models.
- 🔍 **Hyperparameter Tuning** — Use `GridSearchCV` to investigate Ridge regularization strengths.
- 🛡️ **Regularization Analysis** — Compare ordinary least squares with L2 and combined L1/L2 regularization.
- 📏 **Performance Evaluation** — Assess model performance using MAE, RMSE, and R².
- 📈 **Prediction Analysis** — Compare actual and predicted rental counts and identify areas for improvement.

---

## 🧭 Methodology

The notebook follows a structured workflow from data inspection to final evaluation.

| **Step** | **Stage** | **Description** |
| --- | --- | --- |
| **01** | 🚲 Dataset Loading | Load the hourly UCI Bike Sharing data and inspect its structure. |
| **02** | 🔎 Exploratory Data Analysis | Examine the target distribution, hourly rental patterns, and the relationship between temperature and rentals. |
| **03** | 🧹 Feature and Target Preparation | Use `cnt` as the target; exclude `casual` and `registered` because they directly contribute to `cnt`, and exclude `instant` as a record index. |
| **04** | 🕒 Chronological Split | Sort observations by date and use the first 80% for training and the final 20% for testing. |
| **05** | 🧩 Missing-Value Handling | Set up mean imputation for numerical features and mode imputation for categorical features using training data. |
| **06** | ⚙️ Feature Engineering | One-hot encode categorical variables and standardize numerical features. |
| **07** | 📊 Baseline Modeling | Train Linear Regression, Ridge Regression, and Elastic Net. |
| **08** | 🔍 Ridge Tuning | Compare Ridge `alpha` values with `GridSearchCV` and time-series cross-validation. |
| **09** | 📈 Final Evaluation | Evaluate the tuned Ridge model on the held-out test set and visualize actual versus predicted rentals. |

---

## 📊 Regression Models

The project compares three regression approaches to model hourly bike rental counts.

### 🔑 Model Components

| **Model / Tool** | **Function** |
| --- | --- |
| 📈 **Linear Regression** | Establishes a baseline by fitting a linear relationship between input features and rental demand. |
| 🛡️ **Ridge Regression** | Applies L2 regularization to penalize large coefficients and help control model complexity. |
| ⚖️ **Elastic Net Regression** | Combines L1 and L2 regularization. |
| 🔍 **GridSearchCV** | Searches over Ridge regularization strengths using cross-validation. |
| 🕒 **TimeSeriesSplit** | Preserves temporal ordering when forming cross-validation folds. |
| 🧩 **SimpleImputer** | Provides mean imputation for numerical features and mode imputation for categorical features. |
| 🔠 **One-Hot Encoding** | Converts categorical features into model-compatible indicator variables. |
| 📏 **StandardScaler** | Standardizes numerical features before regression modeling. |

> **Modeling strategy:** Establish baseline performance, compare regularized models, tune Ridge Regression, and evaluate the selected model on a chronological hold-out test set.

**Preprocessing note:** The notebook prepares preprocessing matrices before Ridge cross-validation. As documented in the notebook, a more rigorous leakage-safe workflow would put preprocessing and the estimator together in a scikit-learn `Pipeline` inside cross-validation.

---

## 📊 Evaluation Metrics

| **Metric** | **Purpose** |
| --- | --- |
| 📏 **Mean Absolute Error (MAE)** | Average absolute prediction error, measured in rentals. |
| 📉 **Root Mean Squared Error (RMSE)** | Square root of the average squared prediction error; more sensitive to large errors than MAE. |
| 📈 **Coefficient of Determination (R²)** | Measures the proportion of target variation explained by the model relative to a mean-prediction baseline. |

Lower MAE and RMSE indicate smaller errors; higher R² indicates that the model explains more of the variation in the target.

### Baseline Model Results

Held-out test-set results reported in the notebook:

| **Model** | **MAE** | **RMSE** | **R²** |
| --- | ---: | ---: | ---: |
| Linear Regression | 98.798 | 133.835 | 0.632 |
| Ridge Regression | 98.875 | 133.910 | 0.631 |
| Elastic Net | 121.827 | 163.621 | 0.449 |

### Tuned Ridge Results

| **Metric** | **Test Result** |
| --- | ---: |
| MAE | 98.813 rentals |
| RMSE | 133.839 rentals |
| R² | 0.632 |

**Result interpretation:** Linear Regression and Ridge Regression performed very similarly on the held-out test set, explaining approximately **63.2%** of the variation in hourly rental counts. The tuned Ridge model did not meaningfully improve on the baseline. The notebook’s actual-versus-predicted plot also indicates that peak-demand periods remain challenging to predict.

---

## 🛠️ Technologies & Tools

| **Technology** | **Application** |
| --- | --- |
| 🐍 **Python** | Core programming language |
| 🐼 **Pandas** | Data manipulation and inspection |
| 🔢 **NumPy** | Numerical computations |
| 🤖 **Scikit-learn** | Regression models, preprocessing, cross-validation, hyperparameter tuning, and evaluation |
| 📈 **Matplotlib** | Data and prediction visualizations |
| 📊 **Seaborn** | Exploratory data analysis and statistical visualization |
| 📓 **Jupyter Notebook** | Interactive analysis and documentation |

---

## 📁 Project Structure

```text
Bike_Rental_Regression/
├── README.md
├── Bike_Re_submission_ready.ipynb
├── Bike_Re_submission_ready.html
├── requirements.txt          # Optional: Python dependencies
├── figures/                  # Optional: exported plots
└── results/                  # Optional: saved evaluation results
```

- **`README.md`** — Project overview, methodology, tools, and results.
- **`Bike_Re_submission_ready.ipynb`** — Cleaned notebook containing the analysis, models, evaluation, and visualizations.
- **`Bike_Re_submission_ready.html`** — Browser-friendly export of the notebook.
- **`requirements.txt`** — Optional list of Python dependencies for reproducibility.
- **`figures/`** and **`results/`** — Optional folders for exported plots and result files.

The folder structure is illustrative; adjust filenames to match the files actually committed to your repository. The notebook loads the UCI Bike Sharing hourly dataset, so ensure the dataset is available at the path expected by the notebook.

---

## 🚀 Key Findings

- ✅ The original dataset contains **17,379 hourly observations** and **17 columns**, with no missing values reported in the notebook.
- 🕒 The data is evaluated chronologically, with the final 20% held out for testing.
- 📈 Rental activity is lower overnight and tends to rise during commuting hours.
- 🌡️ Warmer temperatures are generally associated with more rentals, though demand depends on multiple factors.
- 📊 Linear Regression and Ridge Regression achieved similar test performance, with R² around **0.63**.
- ⚖️ Elastic Net performed worse than the other two baseline models on this split.
- 🔍 Ridge tuning did not produce a meaningful improvement on the held-out test set.
- 📌 Potential next steps include evaluating nonlinear models, improving feature engineering, and moving all preprocessing into a pipeline used within cross-validation.
- ⚠️ Regression predictions can be negative even though rental counts cannot; inspect or appropriately constrain predictions before operational use.

---

## 🔬 Scientific Motivation

Bike rental demand varies with time, weather, and other contextual factors. Understanding these relationships can support demand planning and help identify periods when rental activity is likely to be high or low.

**Regression provides a data-driven approach to estimating hourly rental demand from observed features.** Comparing ordinary Linear Regression with Ridge and Elastic Net illustrates how regularization affects performance on structured data. Chronological splitting is especially important because future observations should not be used to train a model that is evaluated on earlier observations.

This project demonstrates a practical regression workflow while highlighting the importance of avoiding target leakage, respecting temporal order, evaluating on held-out data, and recognizing model limitations.

---

## 📌 Project Focus

**Bike Sharing Data → Exploratory Analysis → Chronological Split → Preprocessing → Linear Regression → Ridge and Elastic Net → Cross-Validation → Test Evaluation → Prediction Analysis**
