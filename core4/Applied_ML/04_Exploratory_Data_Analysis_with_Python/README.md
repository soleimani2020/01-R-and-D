# Project 04: Exploratory Data Analysis and Statistics with Python

> 📊 **Exploratory Data Analysis** · 📐 **Statistical Analysis** · 📈 **Data Visualization** · 🐍 **Python**

---

## 📌 Project Overview

This project explores the application of **Python for exploratory data analysis (EDA) and statistical investigation**. It focuses on understanding dataset structure, assessing data characteristics, identifying patterns, and examining relationships between variables.

The project combines probability simulation and the analysis of a YouTube dataset containing video metadata and engagement metrics. The workflow includes data preparation, descriptive statistics, visualization, correlation analysis, and hypothesis testing to extract meaningful insights from data.

Rather than immediately building predictive models, the emphasis is on understanding the data, formulating relevant questions, and drawing conclusions supported by observations and appropriate statistical evidence.

### 🎯 Main Goal

Develop a systematic and reproducible workflow for exploring datasets, applying statistical methods, visualizing patterns, and communicating reliable, data-supported insights.

---

## 🎯 Project Objectives

- 🎲 **Probability Simulation** — Generate random samples from a normal distribution and investigate their empirical properties.
- 📂 **Data Exploration** — Inspect dataset dimensions, feature types, unique values, and overall structure.
- 🧹 **Data Quality Assessment** — Examine data types and identify potential inconsistencies.
- 📊 **Descriptive Statistics** — Summarize numerical variables, particularly video engagement metrics.
- 📈 **Distribution Analysis** — Investigate variable distributions, variability, skewness, and potential deviations from normality.
- 🔗 **Relationship Analysis** — Explore associations between engagement metrics and video categories.
- 🧪 **Statistical Investigation** — Apply suitable methods, including chi-square and Mann–Whitney U tests, when justified by the research question.
- 💡 **Insight Generation** — Interpret findings, acknowledge limitations, and communicate evidence-based conclusions.

---

## 🧭 Methodology

The project follows a structured analytical workflow, progressing from probability simulation and initial data inspection to statistical interpretation.

| Step | Stage | Description |
|:---:|---|---|
| **01** | 🎲 Probability Simulation | Generate samples from a normal distribution and compare empirical statistics with theoretical parameters. |
| **02** | 📂 Data Loading | Import the YouTube dataset and inspect its dimensions, columns, and data types. |
| **03** | 🔎 Initial Exploration | Examine category frequencies and summarize the dataset's characteristics. |
| **04** | 🧹 Data Preparation | Convert publication dates to datetime format and filter selected video categories. |
| **05** | 📊 Descriptive Statistics | Calculate summary measures for engagement variables such as likes and views. |
| **06** | 📈 Temporal Analysis | Investigate monthly changes in average likes and video publication counts across categories. |
| **07** | 🔗 Relationship Analysis | Calculate correlations between views, likes, dislikes, and comments. |
| **08** | 🧪 Statistical Investigation | Apply a chi-square test of independence and a Mann–Whitney U test to investigate selected relationships. |
| **09** | 💡 Interpretation and Reporting | Summarize findings, discuss limitations, and communicate data-supported conclusions. |

---

## 📊 Statistical Methods

The project explores fundamental statistical techniques for summarizing data, understanding variability, and investigating relationships.

| Method | Application |
|---|---|
| 📏 **Measures of Central Tendency** | Describe typical values using the mean and median. |
| 📐 **Measures of Dispersion** | Quantify variability using standard deviation and quantiles. |
| 📊 **Distribution Analysis** | Investigate the shape, spread, and skewness of numerical variables such as view counts. |
| 🔗 **Correlation Analysis** | Examine the strength and direction of associations between views, likes, dislikes, and comments. |
| 🧪 **Chi-Square Test of Independence** | Investigate whether video category and publication month are statistically associated. |
| 📉 **Mann–Whitney U Test** | Compare the distributions of likes for videos with comments disabled and enabled. |
| 📈 **Time-Based Grouping** | Examine monthly trends in engagement metrics and video publication counts across categories. |

> **Important distinction:** Statistical association does not necessarily imply causation. A significant test result should be interpreted in the context of the dataset, study design, and assumptions of the method.

---

## 🛠️ Technologies & Tools

| Technology | Application |
|---|---|
| 🐍 **Python** | Core programming language for data analysis |
| 🐼 **Pandas** | Data manipulation, cleaning, grouping, and exploratory analysis |
| 🔢 **NumPy** | Numerical computations and probability simulation |
| 📈 **Matplotlib** | Visualization of distributions, trends, and relationships |
| 🎨 **Seaborn** | Statistical visualization and exploration of variable relationships |
| 📐 **Statsmodels** | Statistical graphics, including Q–Q plots |
| 🧪 **SciPy** | Statistical hypothesis testing |
| 📓 **Jupyter Notebook** | Interactive analysis, experimentation, and documentation |

---

## 📁 Project Structure

```text
04_EDA_and_Statistics/
├── README.md
├── project.ipynb
└── youtube.csv
```

- **`README.md`** — Project overview, objectives, methodology, tools, and documentation.
- **`project.ipynb`** — Interactive notebook containing probability simulation, data exploration, visualizations, statistical analysis, and findings.
- **`youtube.csv`** — YouTube video metadata and engagement data used in the analysis, if permitted to include it in the repository.

The dataset does not necessarily need to be included in the repository. If redistribution is not permitted, document its source and explain how to obtain it.

---

## 🚀 Expected Outcomes

The project aims to establish a reproducible and systematic approach to exploratory data analysis and statistical investigation.

Key outcomes include:

- ✅ An empirical investigation of normal-distribution sampling.
- 📊 Descriptive summaries of YouTube engagement metrics.
- 📈 Visualizations of monthly trends across selected video categories.
- 🔗 Correlation analysis of views, likes, dislikes, and comments.
- 📉 An assessment of whether view counts follow a normal distribution.
- 🧪 Statistical tests investigating associations involving category, publication month, and comment settings.
- 💡 Clear interpretation of findings, including analytical limitations.
- 📌 A documented foundation for subsequent statistical modeling or machine learning.

Actual findings should be reported using the notebook's calculated results rather than assumed in advance.

---

## 🔬 Scientific Motivation

Exploratory data analysis is a fundamental stage of data science. Before developing predictive models or drawing substantive conclusions, analysts need to understand the structure, quality, variability, and limitations of the available data.

**Descriptive statistics summarize observed patterns, while statistical methods help evaluate whether those patterns provide meaningful evidence for a particular question.** Probability simulation provides an additional way to explore the relationship between theoretical distributions and empirical samples.

This project develops practical skills in statistical reasoning, data visualization, hypothesis testing, and evidence-based interpretation using real-world video engagement data.

---

## 📌 Project Focus

**Probability Simulation → Data Exploration → Data Preparation → Descriptive Statistics → Temporal Analysis → Correlation Analysis → Hypothesis Testing → Insight Generation**
