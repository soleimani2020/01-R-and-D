AdaBoost From Scratch in Python 🚀
Overview

This repository presents a from-scratch implementation of AdaBoost (Adaptive Boosting) using Python and NumPy, without relying on external machine learning libraries such as scikit-learn.

AdaBoost is one of the foundational ensemble learning algorithms that transforms multiple weak learners into a powerful classifier by iteratively focusing on previously misclassified samples.

Core Concepts Covered:
Weak learner optimization (decision stumps)
Weighted sample training
Error minimization
Alpha (learner importance) calculation
Iterative boosting
Ensemble classifier construction
Key Features ✨
📌 Pure Python + NumPy implementation
🌲 Decision stump weak learners
🎯 Threshold and polarity optimization
⚖️ Weighted error minimization
📈 Adaptive sample weight updates
🔁 Multiple boosting iterations
🧠 Educational and interpretable design
How AdaBoost Works ⚙️

For each boosting round:

1️⃣ Initialize Sample Weights

All samples begin with equal importance:

w
i
	​

=
N
1
	​


Where:

N = number of training samples
2️⃣ Train a Weak Learner

Search for the best decision stump by optimizing:

Feature
Threshold
Polarity
3️⃣ Compute Weighted Error
error=
i=1
∑
N
	​

w
i
	​

⋅1(y
i
	​


=h
i
	​

(x
i
	​

))

Where:

w_i = sample weight
y_i = true label
h_i(x_i) = weak learner prediction
\mathbb{1} = indicator function
4️⃣ Compute Learner Importance
α=
2
1
	​

ln(
error
1−error
	​

)
Interpretation:
Lower error → higher alpha
Better classifiers receive greater influence
5️⃣ Update Sample Weights
w
i
	​

←w
i
	​

exp(−αy
i
	​

h
i
	​

(x
i
	​

))
6️⃣ Normalize Weights
w
i
	​

←
∑
j=1
N
	​

w
j
	​

w
i
	​

	​

7️⃣ Repeat

Repeat the process for:

n_clf

boosting rounds.

Final Strong Classifier 🏆

The final prediction is:

H(x)=sign(
t=1
∑
T
	​

α
t
	​

h
t
	​

(x))

Where:

T = number of weak learners
\alpha_t = learner weight
h_t(x) = learner prediction
Example Dataset 📊
import numpy as np

X = np.array([
    [1, 2],
    [2, 3],
    [3, 4],
    [4, 5]
])

y = np.array([1, 1, -1, -1])
Example Usage 🖥️
clfs = adaboost_fit(X, y, n_clf=3)
print(clfs)
Example Output 📌
[
    {'feature': 0, 'threshold': 3, 'polarity': 1, 'alpha': ...},
    {'feature': 1, 'threshold': 4, 'polarity': -1, 'alpha': ...}
]
Repository Structure 📂
P36/
├── adaboost.py
└── README.md
Installation 🔧
git clone https://github.com/soleimani2020/01-R-and-D.git
cd "01-R-and-D/core3/Deep_ML/Machine Learning/P36"
pip install numpy
Usage ▶️

Run the implementation:

python adaboost.py
Learning Objectives 🎓

This project strengthens understanding of:

Ensemble learning
Boosting algorithms
Decision stumps
Weighted optimization
Machine learning mathematics
Classifier theory
Future Enhancements 🔮
✅ Prediction on unseen data
📉 Decision boundary visualization
📊 Accuracy and precision metrics
⚔️ Comparison with scikit-learn
🌍 Real-world datasets
⚡ Performance optimization
Why This Project Matters 🌟

Building AdaBoost from scratch provides deeper insight into:

How boosting improves weak learners
The mathematical intuition behind adaptive weighting
Ensemble classifier design
Core machine learning optimization principles

This repository is ideal for:

Students learning ML fundamentals
Researchers exploring boosting algorithms
Developers interested in interpretable ML systems
Author 👨‍🔬

Alireza Soleimani
PhD in Biophysics | Machine Learning | Algorithm Development

GitHub: https://github.com/soleimani2020
License 📜

This project is licensed under the MIT License.

Connect 🌐

If you found this project useful:

⭐ Star the repository
🔁 Share with others
💡 Contribute improvements
