# Project 02: Bicycle vs. Motorcycle Image Classification with CNNs

> 🖼️ **Computer Vision** · 🧠 **Deep Learning** · 🔍 **Image Recognition** · 🤖 **Convolutional Neural Networks**

---

## 📌 Project Overview

Image classification is a fundamental task in computer vision that involves assigning predefined categories to images based on their visual content.

This project explores the development of a deep learning model using **Convolutional Neural Networks (CNNs)** to distinguish between images of bicycles and motorcycles. CNNs are particularly effective for image analysis because they learn spatial features such as edges, textures, shapes, and increasingly complex visual patterns directly from image data.

The project uses TensorFlow and Keras to build, train, and evaluate a CNN-based image classification model. The workflow includes image preprocessing, model architecture design, training optimization, performance evaluation, and visualization of predictions.

The project also explores techniques intended to improve generalization and reduce overfitting, including batch normalization, dropout, data augmentation, and early stopping.

---

## 🎯 Main Goal

Develop and evaluate a CNN-based model capable of distinguishing between **bicycles and motorcycles** by learning meaningful visual representations from image data.

---

## 🎯 Project Objectives

- 🖼️ **Image Preprocessing:** Prepare image data through normalization and appropriate formatting.
- 🔢 **Label Encoding:** Represent the two target classes in a format compatible with the selected loss function.
- 🧠 **CNN Architecture:** Build a neural network using convolutional and pooling layers.
- 🔍 **Feature Extraction:** Learn spatial representations from image pixels.
- ⚙️ **Model Optimization:** Train the network using an appropriate optimizer and loss function.
- 🛡️ **Regularization:** Investigate dropout, batch normalization, and data augmentation.
- 📉 **Training Monitoring:** Track training and validation loss and accuracy.
- 📏 **Performance Evaluation:** Evaluate classification performance on unseen images.
- 📊 **Visualization:** Analyze training curves, predictions, and classification errors.
- 🔬 **Model Interpretation:** Explore techniques such as Grad-CAM or occlusion sensitivity when appropriate.

---

## 🧭 Methodology

| Step | Description | Purpose |
|---|---|---|
| 1. Data Loading | Load CIFAR-100 images and select the bicycle and motorcycle classes | Establish the training and evaluation data |
| 2. Preprocessing | Normalize pixel values and prepare binary class labels | Make the data suitable for neural network training |
| 3. Data Augmentation | Apply suitable transformations to training images | Improve robustness and generalization |
| 4. CNN Architecture | Define convolutional, activation, pooling, and dense layers | Learn increasingly complex visual features |
| 5. Model Compilation | Select an optimizer, loss function, and evaluation metrics | Configure the training process |
| 6. Model Training | Train the network using mini-batches over multiple epochs | Learn image-to-class mappings |
| 7. Regularization | Apply techniques such as dropout and batch normalization | Help control overfitting |
| 8. Training Monitoring | Track loss and accuracy across training and validation data | Assess convergence and generalization |
| 9. Model Evaluation | Evaluate performance on held-out test images | Measure predictive performance |
| 10. Error Analysis | Examine incorrect predictions between the two classes | Identify model limitations |

---

## 🧠 CNN Architecture

A Convolutional Neural Network learns spatial features through successive layers.

| Layer | Function |
|---|---|
| Convolutional Layer | Applies learnable filters to detect local visual patterns |
| Activation Function | Introduces nonlinearity, commonly using ReLU |
| Max Pooling | Reduces spatial dimensions while retaining prominent activations |
| Batch Normalization | Normalizes intermediate activations during training |
| Dropout | Randomly deactivates selected units during training to help reduce overfitting |
| Flatten or Global Pooling | Converts spatial feature maps into a representation suitable for classification |
| Dense Layer | Combines learned features to produce class predictions |
| Sigmoid Output | Produces a probability for the binary classification task |

The final architecture should be adapted to the dataset size, image resolution, and available computational resources.

---

## 🗂️ Dataset

The project uses **CIFAR-100**, a labeled image dataset containing 100 classes of small color images measuring **32 × 32 pixels**.

This implementation selects two classes:

- **Bicycle**
- **Motorcycle**

CIFAR-100 contains 600 images per class, with 500 training images and 100 test images per class in its standard split. The selected binary dataset therefore contains 1,000 training images and 200 test images before creating a validation subset.

Potential preprocessing steps include:

- Normalizing pixel values to a suitable numerical range.
- Mapping the bicycle and motorcycle labels to binary targets.
- Creating a validation set from the training data.
- Applying augmentation to training images where appropriate.
- Keeping the test set separate from model training and hyperparameter tuning.

The dataset's source and applicable license should be documented in the project README.

---

## 🛡️ Regularization and Generalization

Deep neural networks can overfit their training data, particularly when model capacity is large relative to the available training examples.

Potential strategies include:

- **Data Augmentation:** Generate transformed training images to expose the model to additional visual variations.
- **Dropout:** Reduce reliance on specific neural network units.
- **Batch Normalization:** Help stabilize intermediate activations during training.
- **Early Stopping:** Stop training when validation performance no longer improves according to the chosen criterion.
- **Architecture Selection:** Balance model complexity with dataset size and computational cost.

These techniques can improve generalization, but their effectiveness should be verified experimentally.

---

## 📏 Model Evaluation

Classification performance can be evaluated using several complementary metrics.

| Metric | Description |
|---|---|
| Accuracy | Proportion of images assigned the correct class |
| Training Loss | Measures the model's objective on training data |
| Validation Loss | Helps monitor generalization during model development |
| Precision | Measures the proportion of predicted positive-class assignments that are correct |
| Recall | Measures the proportion of actual positive-class instances correctly identified |
| F1-Score | Combines precision and recall |
| Confusion Matrix | Shows how often bicycles and motorcycles are correctly or incorrectly classified |

For this binary classification task, the positive class should be defined clearly—for example, **motorcycle = 1** and **bicycle = 0**—so that precision and recall can be interpreted consistently.

The final performance assessment should use a held-out test set that was not used for model selection.

---

## 🔬 Model Interpretation

Beyond overall accuracy, visualization can help investigate what the model learns and where it makes mistakes.

Potential techniques include:

- **Prediction Visualization:** Display images alongside predicted and actual labels.
- **Confusion Matrix:** Identify which class is more frequently misclassified.
- **Grad-CAM:** Highlight image regions that contribute to a selected class prediction, when compatible with the architecture.
- **Occlusion Sensitivity:** Examine how predictions change when image regions are masked.

These methods can provide useful diagnostic evidence, but their outputs should not be interpreted as definitive explanations of the model's internal reasoning.

---

## 🛠️ Technologies & Tools

| Category | Technologies |
|---|---|
| Programming Language | Python |
| Deep Learning Framework | TensorFlow, Keras |
| Numerical Computing | NumPy |
| Data Analysis | Pandas |
| Visualization | Matplotlib, Seaborn |
| Model Evaluation | Scikit-learn |
| Development Environment | Jupyter Notebook |

---

## 📁 Project Structure

```text
02_CNN_Bicycle_vs_Motorcycle/
│
├── README.md
├── project.ipynb
│
├── figures/               # Training curves and evaluation plots
├── models/                # Saved model files, if included
└── results/               # Evaluation outputs, if included
```

The folder structure can be adapted to the final implementation. Large datasets and model files do not necessarily need to be included in the repository.

---

## 🚀 Expected Outcomes

The project aims to produce:

- A working CNN architecture for binary image classification.
- A reproducible image preprocessing and training workflow.
- Training and validation curves for monitoring model performance.
- Quantitative evaluation on held-out test data.
- Visualization of predictions and classification errors.
- An assessment of overfitting, generalization, and potential improvements.

Actual performance should be reported from the trained model rather than assumed in advance.

---

## 🔬 Scientific Motivation

CNNs have played an important role in modern computer vision by enabling neural networks to learn hierarchical spatial representations directly from images.

Their ability to combine local feature extraction with nonlinear transformations makes them useful for recognizing visual differences between objects such as bicycles and motorcycles.

This project provides practical experience with deep learning, model optimization, regularization, and evaluation while illustrating the challenges of building reliable image classification systems.

---

## 📌 Project Focus

**CIFAR-100 → Class Selection → Preprocessing → CNN Architecture → Training and Optimization → Evaluation → Error Analysis → Model Interpretation**
