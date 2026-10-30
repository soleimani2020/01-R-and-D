# Who Will Repay a Loan? Predicting Risk and Making the Approval Decision

> 💳 **Credit Risk Analysis** · 🎯 **Classification** · 💰 **Cost-Sensitive Decisions** · 🤖 **Machine Learning**

---

## 📌 Project Overview

Credit risk screening is an important application of machine learning in financial services. A lender must decide, applicant by applicant, whether to approve a loan, and the two possible mistakes are not equally costly. Approving a borrower who later defaults can lose most of the principal, while turning away a borrower who would have repaid only forgoes the interest.

This project uses the public German Credit dataset (1,000 past applicants, 20 recorded attributes, outcome good/bad) to build a risk score, check that its probabilities can be trusted, and turn it into an approve/reject rule that reflects the real cost of each mistake.

The workflow covers exploratory data analysis, preprocessing, model comparison with repeated cross-validation, calibration analysis, cost-based cut-off selection, a final check on held-out applicants, and a review of feature importance and age-group outcomes.

The objective is to show that a good ranking model is not enough: the approval cut-off must reflect the unequal costs of errors.

---

## 🎯 Main Goal

Develop and evaluate a loan screening model that ranks applicants by default risk, produces trustworthy probabilities, and supports an approval rule that minimizes expected cost.

---

## 🎯 Project Objectives

- 📊 **Exploratory Data Analysis:** Investigate borrower characteristics, loan attributes, and repayment outcomes.
- 🧹 **Data Inspection:** Check data types, missing values, and the target distribution.
- ⚙️ **Feature Engineering:** Add a `monthly_payment` feature (loan amount divided by duration) as an affordability signal.
- ⚖️ **Class Imbalance Analysis:** Account for the 30% default rate with stratified splitting and stratified cross-validation.
- 🎯 **Loan Classification:** Train models to predict the probability of default.
- 🔍 **Model Comparison:** Evaluate four candidates under identical preprocessing and validation.
- 📏 **Performance Evaluation:** Assess ranking quality (ROC AUC, average precision) and probability quality (Brier score).
- 📈 **Calibration Analysis:** Check whether a predicted 10% risk really means about 10 in 100 defaults.
- 💰 **Cost-Based Decision Rule:** Find the cut-off that minimizes total cost and compare it with the naive 50% rule.
- 🔎 **Model Interpretation:** Use permutation importance to find which details drive the score.
- ⚖️ **Responsible Modeling:** Check outcomes across age groups and discuss fairness and limitations.

---

## 🧭 Methodology

| Step | Description | Purpose |
|---|---|---|
| 1. Data Collection | Load `german_credit.csv`, or fetch the dataset from OpenML | Establish the analytical dataset |
| 2. Data Inspection | Examine data types, missing values, and target distribution | Assess data quality |
| 3. Exploratory Analysis | Default rate by checking-account status, loan size and duration by outcome | Understand what separates good from bad outcomes |
| 4. Feature Engineering | Add `monthly_payment` | Improve the representation of affordability |
| 5. Data Splitting | Stratified 75/25 train/test split | Evaluate generalization, keep the test set untouched |
| 6. Baseline Modeling | Do-nothing classifier predicting the class prior | Provide a reference floor |
| 7. Model Training | Logistic regression, random forest, gradient boosting | Learn patterns associated with default |
| 8. Model Evaluation | Stratified 5-fold cross-validation repeated 3 times | Get a stable estimate on a small dataset |
| 9. Calibration Check | Out-of-fold predictions and calibration curves | Test whether probabilities can be acted on |
| 10. Cut-off Selection | Try cut-offs from 5% to 95% on out-of-fold training predictions | Find the lowest-cost rule |
| 11. Final Check | Score the 250 held-out applicants once, compare four policies | Confirm the rule saves money on unseen data |
| 12. Interpretation | Permutation importance and age-group analysis | Understand drivers and check for large disparities |

---

## 📊 Exploratory Data Analysis

- **Target Distribution:** 30% of applicants defaulted.
- **Checking-Account Status:** Default rate rises sharply with account status, from about 12% (no checking account) to about 49% (overdrawn account).
- **Loan Duration:** Defaulters typically had longer loans (median about 24 months vs 18).
- **Loan Amount:** Defaulters had somewhat larger loans in the upper range, though the medians are similar.
- **Missing Values:** None.

These analyses reveal associations in the observed data, but they do not show that a characteristic causes default.

---

## ⚙️ Feature Engineering and Preprocessing

- **Feature Scaling:** Numerical features are standardized.
- **Categorical Encoding:** Categorical features are one-hot encoded.
- **Feature Construction:** `monthly_payment` is added.
- **Pipeline Construction:** Preprocessing and model sit in one Scikit-learn pipeline, so everything that learns from data is re-fitted on each training fold and never sees held-out data.

---

## 🎯 Classification Models

| Algorithm | Role in this project |
|---|---|
| Baseline (class prior) | Floor that every other model must beat |
| Logistic Regression | Interpretable model, **selected as the final model** |
| Random Forest | Captures nonlinear interactions |
| Gradient Boosting | Flexible ensemble |

Model selection used ranking quality, calibration around the cut-off, and simplicity of explanation.

---

## ⚖️ Class Imbalance and Cost-Sensitive Decisions

The two errors have different costs: approving a borrower who defaults costs **5 units**, while rejecting a borrower who would have repaid costs **1 unit**. This 5-to-1 ratio is the one published with the German Credit dataset.

Approving is worth it only when 5 × P(default) < 1 × (1 − P(default)), which gives a cut-off of 1/6 (about 17%). The notebook compares this theoretical cut-off with the best cut-off found on the training data and with the naive 50% rule.

The cut-off is chosen on training data only. The test set is used once, at the end.

---

## 📏 Model Evaluation

| Metric | Why it is used |
|---|---|
| ROC AUC | Does a random defaulter get a higher score than a random repayer? |
| Average Precision | How concentrated are defaulters at the top of the ranking? |
| Brier Score | Quality of the predicted probabilities (lower is better) |
| Calibration Curve | Does predicted risk match observed default frequency? |
| Total Cost | Cost of both mistakes combined, which is the decision metric |

Cross-validation results (training data, 15 test scores per model):

| Model | ROC AUC | Avg Precision | Brier |
|---|---|---|---|
| Logistic regression | 0.774 | 0.589 | 0.173 |
| Random forest | 0.768 | 0.600 | 0.171 |
| Gradient boosting | 0.752 | 0.574 | 0.190 |
| Baseline | 0.500 | 0.300 | 0.210 |

The three real models are effectively tied on ranking (the fold-to-fold spread is about 0.03). Calibration decided the choice: logistic regression stayed closest to the diagonal around the 17% cut-off, gradient boosting was too optimistic at low predicted risk, and random forest was erratic.

---

## 💰 Final Result on Held-Out Applicants

Total cost on the 250 test applicants:

| Policy | Total cost | Per applicant |
|---|---|---|
| **Model, chosen cut-off (1/6)** | **131** | **0.52** |
| Reject everyone | 175 | 0.70 |
| Model, naive 50% rule | 190 | 0.76 |
| Approve everyone | 375 | 1.50 |

The cost-aware rule is about 31% cheaper than the naive 50% rule. The naive rule is even worse than rejecting everyone, so a good ranking model only pays off when its cut-off reflects the real costs. The best training cut-off (0.20) performed the same as the theoretical 1/6 (130 vs 131), so the simple rule "reject above 1/6" is enough.

---

## 🔎 Model Interpretation and Risk Analysis

- **Permutation Importance:** The model relies most on checking-account status (AUC drop 0.094), then loan duration (0.045), purpose (0.034), savings status (0.030), and credit history (0.029). Loan amount matters less (0.013).
- **Age Groups:** The rule rejects about 60% of applicants against a 30% default rate, because a cost-aware cut-off is deliberately cautious. Reject rates follow real risk: highest under 25 (72% rejected, 40% default) and lowest at ages 35-49 (49% rejected, 19% default). The gap between reject rate and default rate is about 30 percentage points in the three largest bands. The 50+ band has only 22 applicants and is too noisy to interpret.

Feature importance describes model behavior; it does not establish causality. Age is a model input here, which would raise legal and fairness questions in real lending.

---

## 🛠️ Technologies & Tools

| Category | Technologies |
|---|---|
| Programming Language | Python 3.10+ |
| Data Manipulation | Pandas, NumPy |
| Data Visualization | Matplotlib (3.9+) |
| Machine Learning | Scikit-learn |
| Development Environment | Jupyter Notebook |

```bash
pip install numpy pandas matplotlib scikit-learn jupyter
```

---

## 📁 Project Structure

```text
credit-risk-screening/
│
├── README.md
├── credit_risk_screening.ipynb
└── german_credit.csv          # Optional local copy of the data
```

If `german_credit.csv` is not present, the notebook downloads the dataset from OpenML (internet connection required).

---

## 🚀 How to Run

```bash
git clone <your-repo-url>
cd <your-repo-name>
jupyter notebook credit_risk_screening.ipynb
```

Use **Run All**. To test a different cost ratio, change `COST_BAD_APPROVED` and `COST_GOOD_REJECTED` in the setup cell, and update `chosen_cutoff` in the final policy cell to cost_good_rejected / (cost_bad_approved + cost_good_rejected).

---

## ⚠️ Limitations

- **Small data:** With 1,000 applicants (250 in the test set), the exact savings are approximate. The order of the policies is clear, but the percentages could shift with a different split.
- **Assumed costs:** The 5-to-1 ratio comes with the dataset. A real lender should use its own figures for losses and interest.
- **Old, selected data:** The data are historical German loans, and outcomes are only known for approved applicants, so the model may not generalize to people the bank would have rejected.
- **Cost of caution:** The rule rejects many applicants who would have repaid. Its cost counts only the forgone interest.
- **Fairness:** Age is an input to the model, and the notebook only checks outcomes across age groups. Using age in real lending decisions is restricted or illegal in many countries.

---

## 🔬 Scientific Motivation

Credit risk analysis is an important application of statistical modeling and machine learning. However, financial risk prediction involves more than maximizing ranking accuracy. Unequal error costs, probability quality, small samples, and fairness considerations all affect whether a model is useful.

This project explores these challenges in a reproducible workflow, connecting predictive modeling with calibration, cost-based decision making, and responsible financial data analysis.

---

## 📚 Data Source

The German Credit dataset is available from the [UCI Machine Learning Repository](https://archive.ics.uci.edu/) and [OpenML](https://www.openml.org/) (data id 31). Check the dataset's terms before redistributing `german_credit.csv`.

---

## 📌 Project Focus

**Loan Data → Exploratory Analysis → Preprocessing → Model Comparison → Calibration → Cost-Based Cut-off → Held-Out Evaluation → Interpretation → Responsible Analysis**
