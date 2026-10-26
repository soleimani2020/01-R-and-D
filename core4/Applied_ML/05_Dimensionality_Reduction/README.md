# Project 05: Dimensionality Reduction for Human Activity Recognition

> 📉 **Dimensionality Reduction** · 📱 **Sensor Data Analysis** · 🤖 **Machine Learning** · 🎨 **Data Visualization**

---

## 📌 Project Overview

This project investigates the application of **dimensionality reduction techniques to Human Activity Recognition (HAR)**. The objective is to explore high-dimensional sensor data and determine how lower-dimensional representations can reveal patterns associated with different human activities.

Using a structured machine learning workflow, the project examines sensor features, prepares the data for analysis, applies dimensionality reduction, and visualizes the resulting representations. Where classification is included, the analysis also investigates whether the reduced feature space retains information relevant to distinguishing between activities.

The project emphasizes the relationship between feature dimensionality, information preservation, data visualization, and downstream predictive performance.

### 🎯 Main Goal

Investigate how dimensionality reduction can transform high-dimensional sensor data into more compact representations while preserving meaningful structure for exploring and distinguishing human activities.

---

## 🎯 Project Objectives

- 📂 **Dataset Exploration** — Examine the dataset structure, sensor features, activity labels, and feature dimensionality.
- 🧹 **Data Preprocessing** — Prepare sensor measurements for analysis and address relevant data quality issues.
- 📊 **Feature Analysis** — Investigate feature scales, distributions, and potential redundancy.
- 📉 **Dimensionality Reduction** — Apply suitable techniques to represent high-dimensional data using fewer dimensions.
- 🎨 **Pattern Visualization** — Explore lower-dimensional representations to identify potential activity-related patterns.
- 🤖 **Classification Analysis** — Evaluate whether reduced representations retain useful information for activity classification, if included.
- ⚖️ **Performance Comparison** — Compare original and reduced feature spaces using appropriate evaluation criteria.
- 🔬 **Result Interpretation** — Assess information retention, computational considerations, and limitations of the selected methods.

---

## 🧭 Methodology

The project follows a structured workflow from data exploration to the interpretation of lower-dimensional representations.

| Step | Stage | Description |
|:---:|---|---|
| **01** | 📂 Dataset Exploration | Load the human activity recognition dataset and inspect its features, labels, and structure. |
| **02** | 🧹 Data Preprocessing | Prepare sensor features and activity labels for analysis. |
| **03** | 📊 Feature Analysis | Examine feature distributions, scales, and potential redundancy. |
| **04** | 📉 Dimensionality Reduction | Transform the original features into a lower-dimensional representation using appropriate techniques. |
| **05** | 🎨 Visualization | Visualize reduced representations and investigate potential patterns across activity classes. |
| **06** | 🤖 Classification Analysis | Evaluate activity classification using original or reduced features, if applicable. |
| **07** | 📈 Performance Evaluation | Compare classification metrics or representation characteristics, depending on the project objectives. |
| **08** | 🔬 Interpretation | Examine information retention, activity-related structure, and methodological limitations. |

---

## 🧠 Dimensionality Reduction

Dimensionality reduction transforms high-dimensional data into a representation with fewer dimensions while attempting to preserve information relevant to the analysis.

### 🔑 Potential Techniques

| Technique | Purpose |
|---|---|
| 📐 **Principal Component Analysis (PCA)** | Projects data onto orthogonal directions that capture the greatest possible variance, enabling compact representations and analysis of explained variance. |
| 🌀 **t-distributed Stochastic Neighbor Embedding (t-SNE)** | Constructs a low-dimensional representation designed to preserve local neighborhood relationships, primarily for visualization. |
| 🧭 **Uniform Manifold Approximation and Projection (UMAP)** | Constructs a lower-dimensional representation intended to preserve aspects of the data's neighborhood structure. |

The techniques serve different purposes and should not be interpreted as interchangeable methods.

- **PCA** provides a linear transformation and supports quantitative analysis of explained variance.
- **t-SNE** is particularly useful for visualizing local relationships, but distances between separated clusters and their apparent sizes require careful interpretation.
- **UMAP** can reveal local structure and may also support transformations of new observations when appropriately configured.

> **Initial strategy:** Establish a baseline using the original feature space, then apply the dimensionality reduction techniques relevant to the project. Evaluate each representation according to its intended purpose rather than assuming that a visually clearer embedding necessarily improves classification.

---

## 📊 Evaluation Strategy

The evaluation depends on whether dimensionality reduction is used primarily for visualization, feature compression, or downstream classification.

| Evaluation Aspect | Purpose |
|---|---|
| 📉 **Dimensionality Reduction** | Determine how many dimensions are retained in the reduced representation. |
| 📐 **Explained Variance** | For PCA, quantify how much of the original variance is captured by the selected components. |
| 🎨 **Visual Structure** | Examine how activity labels are distributed within the reduced feature space. |
| 🤖 **Classification Performance** | Compare suitable classification metrics using original and reduced features, if classification is implemented. |
| ⚖️ **Information Retention** | Investigate whether the reduced representation preserves useful information for the intended task. |

Visual inspection alone does not establish classification quality. When classification is part of the project, quantitative evaluation on appropriately separated data provides stronger evidence of predictive performance.

---

## 🛠️ Technologies & Tools

| Technology | Application |
|---|---|
| 🐍 **Python** | Core programming language |
| 🐼 **Pandas** | Data manipulation and preprocessing |
| 🔢 **NumPy** | Numerical computations and array operations |
| 🤖 **Scikit-learn** | Feature preprocessing, PCA, and classification workflows |
| 🌀 **t-SNE / UMAP** | Nonlinear dimensionality reduction and visualization, where implemented |
| 📈 **Matplotlib** | Visualization of reduced representations and analytical results |
| 🎨 **Seaborn** | Statistical visualization, where appropriate |
| 📓 **Jupyter Notebook** | Interactive analysis, experiments, and documentation |

---

## 📁 Project Structure

```text
05_Human_Activity_Recognition/
├── README.md
└── project.ipynb
```

- **`README.md`** — Project overview, objectives, methodology, technologies, and documentation.
- **`project.ipynb`** — Interactive notebook for dataset exploration, preprocessing, dimensionality reduction, visualization, and evaluation.

---

## 🚀 Expected Outcomes

The project aims to establish a systematic workflow for investigating dimensionality reduction in human activity recognition.

Key outcomes include:

- ✅ A prepared representation of the human activity sensor dataset.
- 📉 Lower-dimensional representations of the original features.
- 🎨 Visual exploration of activity-related patterns and potential class separation.
- 📐 An assessment of information retention using appropriate measures.
- 🤖 A comparison of classification performance, if downstream classification is included.
- 🔍 An understanding of the strengths and limitations of different dimensionality reduction techniques.
- 💡 Evidence-based conclusions about the usefulness of reduced representations for the intended analytical task.

---

## 🔬 Scientific Motivation

Human Activity Recognition commonly relies on sensor measurements that may contain numerous features with different scales, correlations, and levels of redundancy. High-dimensional data can be challenging to visualize and may increase the complexity of subsequent analysis.

**Dimensionality reduction provides a framework for constructing more compact representations while retaining selected characteristics of the original data.** These representations can help reveal activity-related patterns, support exploratory visualization, and investigate whether a smaller feature space is sufficient for downstream machine learning tasks.

This project develops practical skills in feature representation, dimensionality reduction, visualization, and the evaluation of high-dimensional sensor data.

---

## 📌 Project Focus

**Sensor Data → Data Preprocessing → Feature Analysis → Dimensionality Reduction → Visualization → Activity Pattern Analysis → Evaluation → Interpretation**
