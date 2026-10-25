```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.kernel_ridge import KernelRidge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score


RANDOM_STATE = 42
NOISE_LEVEL = 0.05

BLUE = "#0000ff"
RED = "#ff0000"
CURVE_COLORS = [BLUE, RED, "#008000", "#ff00ff", "#000000"]


# Plot settings
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["LM Roman 10", "DejaVu Serif"],
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "mathtext.fontset": "cm",
    "legend.frameon": False,
    "grid.color": "#bfbfbf",
    "grid.linewidth": 0.5,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.04,
    "svg.fonttype": "none",
})


# Normalized form factor of a homogeneous sphere
def sphere_form_factor(q, radius_nm):
    x = np.asarray(q, dtype=float) * (10.0 * radius_nm)

    amplitude = np.ones_like(x)
    mask = np.abs(x) > 1e-4

    z = x[mask]
    amplitude[mask] = 3.0 * (
        np.sin(z) - z * np.cos(z)
    ) / z**3

    z_small = x[~mask]
    amplitude[~mask] = (
        1.0 - z_small**2 / 10.0 + z_small**4 / 280.0
    )

    return amplitude**2


def evaluate(y_true, y_pred):
    return {
        "MAE (nm)": mean_absolute_error(y_true, y_pred),
        "RMSE (nm)": np.sqrt(mean_squared_error(y_true, y_pred)),
        "R²": r2_score(y_true, y_pred),
    }


def add_noise(X, level, seed):
    rng = np.random.default_rng(seed)
    noise = rng.normal(0, level, size=X.shape)
    return np.clip(X * (1 + noise), 0, None)


def axis_limits(values, log=False, padding=0.04):
    values = np.asarray(values, dtype=float).ravel()
    values = values[np.isfinite(values)]

    if log:
        values = np.log10(values[values > 0])

    low, high = values.min(), values.max()
    span = high - low

    if span == 0:
        span = max(abs(low), 1.0) * 0.02

    low -= padding * span
    high += padding * span

    return (10**low, 10**high) if log else (low, high)


def save_figure(fig, filename=None):
    fig.tight_layout(pad=0.5)

    if filename:
        fig.savefig(filename, dpi=300, bbox_inches="tight", pad_inches=0.04)

    plt.show()
    plt.close(fig)


# Generate synthetic scattering data

n_samples = 800
q = np.linspace(0.004, 0.18, 200)

rng = np.random.default_rng(RANDOM_STATE)
y = rng.uniform(3.0, 25.0, n_samples)

X = np.array([sphere_form_factor(q, radius) for radius in y])

print(f"Dataset shape: {X.shape}")
print(f"Radius range: {y.min():.3f}–{y.max():.3f} nm")


# Check the sphere form factor

q_demo = np.linspace(0.001, 0.3, 1000)
P_demo = sphere_form_factor(q_demo, 5.0)

print(f"Expected first minimum: q ≈ {4.493 / 50:.4f} Å⁻¹")

fig, ax = plt.subplots(figsize=(6.4, 3.8))
ax.plot(q_demo, P_demo, color=BLUE, lw=1.4)
ax.set(
    xlabel=r"q (Å$^{-1}$)",
    ylabel="Normalized form factor P(q)",
    title="Sphere form factor (R = 5 nm)",
)
ax.set_yscale("log")
ax.set_xlim(axis_limits(q_demo))
ax.set_ylim(axis_limits(P_demo, log=True))
ax.grid(True)
save_figure(fig, "sans_00_sphere_form_factor.svg")


# Plot representative scattering curves

target_radii = [3, 5, 10, 15, 25]
curve_indices = [
    np.argmin(np.abs(y - radius)) for radius in target_radii
]

fig, ax = plt.subplots(figsize=(6.4, 4.2))

for idx, color in zip(curve_indices, CURVE_COLORS):
    ax.plot(
        q, X[idx], color=color, lw=1.2,
        marker="o", markevery=14, markersize=4,
        markerfacecolor="none",
        label=f"R = {y[idx]:.2f} nm",
    )

ax.set(
    xlabel=r"q (Å$^{-1}$)",
    ylabel="Normalized P(q)",
    title="Noise-free sphere form factors",
)
ax.set_yscale("log")
ax.set_xlim(axis_limits(q))
ax.set_ylim(axis_limits(X[curve_indices], log=True))
ax.grid(True)
ax.legend()
save_figure(fig, "sans_01_noise_free_curves.svg")


# Split the data into training and test sets

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.20, random_state=RANDOM_STATE
)


# Baseline KRR

krr_baseline = KernelRidge(kernel="rbf", alpha=0.01, gamma=0.01)
krr_baseline.fit(X_train, y_train)

y_pred_baseline = krr_baseline.predict(X_test)
baseline_metrics = evaluate(y_test, y_pred_baseline)

print("\nBaseline KRR")
print(pd.Series(baseline_metrics).to_string())


# Plot predicted versus actual radii

limits = axis_limits(np.concatenate([y_test, y_pred_baseline]))

fig, ax = plt.subplots(figsize=(5.2, 5.2))
ax.scatter(
    y_test, y_pred_baseline, s=27,
    facecolors="none", edgecolors=BLUE,
    linewidths=1, label="Test samples",
)
ax.plot(limits, limits, color=RED, linestyle="--", label="Ideal prediction")
ax.set(
    xlabel="Actual radius (nm)",
    ylabel="Predicted radius (nm)",
    title="Baseline KRR",
    xlim=limits,
    ylim=limits,
)
ax.set_aspect("equal", adjustable="box")
ax.grid(True)
ax.legend()
save_figure(fig, "KRR_predicted_vs_actual.svg")


# Tune KRR hyperparameters using five-fold cross-validation

param_grid = {
    "alpha": [1e-4, 1e-3, 1e-2, 1e-1, 1.0],
    "gamma": [1e-4, 1e-3, 1e-2, 1e-1, 1.0],
}

search_clean = GridSearchCV(
    KernelRidge(kernel="rbf"),
    param_grid,
    cv=5,
    scoring="neg_mean_absolute_error",
    n_jobs=-1,
)
search_clean.fit(X_train, y_train)

krr_optimized = search_clean.best_estimator_
y_pred_optimized = krr_optimized.predict(X_test)
optimized_metrics = evaluate(y_test, y_pred_optimized)

print("\nOptimized KRR")
print("Best parameters:", search_clean.best_params_)
print(f"Cross-validation MAE: {-search_clean.best_score_:.6f} nm")

clean_results = pd.DataFrame({
    "Baseline": baseline_metrics,
    "Optimized": optimized_metrics,
})
print(clean_results.to_string(float_format=lambda v: f"{v:.6f}"))


# Compare residuals on clean test data

for predictions, title, filename in [
    (y_pred_baseline, "Baseline KRR", "sans_02_clean_residuals_baseline.svg"),
    (y_pred_optimized, "Optimized KRR", "sans_02_clean_residuals_optimized.svg"),
]:
    residuals = predictions - y_test

    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    ax.scatter(
        y_test, residuals, s=24,
        facecolors="none", edgecolors=BLUE, linewidths=1,
    )
    ax.axhline(0, color=RED, linestyle="--", lw=1.1)
    ax.set(
        xlabel="True radius (nm)",
        ylabel="Predicted − true radius (nm)",
        title=title,
    )
    ax.set_xlim(axis_limits(y_test))
    ax.set_ylim(axis_limits(np.append(residuals, 0)))
    ax.grid(True)
    save_figure(fig, filename)


# Test the clean-trained model on noisy curves

X_test_noisy = add_noise(X_test, NOISE_LEVEL, seed=123)
y_pred_noisy = krr_optimized.predict(X_test_noisy)
noisy_test_metrics = evaluate(y_test, y_pred_noisy)

print("\nOptimized KRR on noisy test data")
print(pd.Series(noisy_test_metrics).to_string())


# Plot noisy scattering curves

test_curve_indices = [
    np.argmin(np.abs(y_test - radius)) for radius in target_radii
]

fig, ax = plt.subplots(figsize=(8.2, 5.2))

for idx, color in zip(test_curve_indices, CURVE_COLORS):
    ax.plot(
        q, X_test_noisy[idx], color=color, lw=1.2,
        label=f"R = {y_test[idx]:.2f} nm",
    )

ax.set(
    xlabel=r"q (Å$^{-1}$)",
    ylabel="Noisy normalized P(q)",
    title="Sphere form factors with 5% relative Gaussian noise",
)
ax.set_yscale("log")
ax.set_xlim(axis_limits(q))
ax.set_ylim(axis_limits(X_test_noisy[test_curve_indices], log=True))
ax.grid(True)
ax.legend()
save_figure(fig, "sans_03_noise_example_combined.svg")


# Train KRR on noisy training curves

X_train_noisy = add_noise(X_train, NOISE_LEVEL, seed=456)

param_grid_noisy = {
    "alpha": [1e-4, 1e-3, 1e-2, 1e-1, 1.0],
    "gamma": [1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0],
}

search_noisy = GridSearchCV(
    KernelRidge(kernel="rbf"),
    param_grid_noisy,
    cv=5,
    scoring="neg_mean_absolute_error",
    n_jobs=-1,
)
search_noisy.fit(X_train_noisy, y_train)

krr_noise_aware = search_noisy.best_estimator_
y_pred_noise_aware = krr_noise_aware.predict(X_test_noisy)
noise_aware_metrics = evaluate(y_test, y_pred_noise_aware)

print("\nNoise-aware KRR")
print("Best parameters:", search_noisy.best_params_)
print(f"Cross-validation MAE: {-search_noisy.best_score_:.6f} nm")


# Compare all models on their respective test conditions

results = pd.DataFrame.from_dict(
    {
        "Baseline | clean → clean": baseline_metrics,
        "Optimized | clean → clean": optimized_metrics,
        "Optimized | clean → noisy": noisy_test_metrics,
        "Noise-aware | noisy → noisy": noise_aware_metrics,
    },
    orient="index",
)

print("\nModel comparison")
print(results.to_string(float_format=lambda v: f"{v:.6f}"))

mae_clean = results.loc["Optimized | clean → noisy", "MAE (nm)"]
mae_aware = results.loc["Noise-aware | noisy → noisy", "MAE (nm)"]

mae_reduction = 100 * (1 - mae_aware / mae_clean)
print(f"\nMAE reduction on noisy test data: {mae_reduction:.2f}%")


# Compare residuals under the same noisy test conditions

residuals_clean = y_pred_noisy - y_test
residuals_aware = y_pred_noise_aware - y_test

for residuals, title, filename in [
    (
        residuals_clean,
        "Model trained on clean curves",
        "sans_04_noise_aware_residuals_clean_trained.svg",
    ),
    (
        residuals_aware,
        "Model trained on noisy curves",
        "sans_04_noise_aware_residuals_noise_aware.svg",
    ),
]:
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    ax.scatter(
        y_test, residuals, s=24,
        facecolors="none", edgecolors=BLUE, linewidths=1,
    )
    ax.axhline(0, color=RED, linestyle="--", lw=1.1)
    ax.set(
        xlabel="True radius (nm)",
        ylabel="Predicted − true radius (nm)",
        title=title,
    )
    ax.set_xlim(axis_limits(y_test))
    ax.set_ylim(axis_limits(np.append(residuals, 0)))
    ax.grid(True)
    save_figure(fig, filename)


# Summarize residual errors

residual_summary = pd.DataFrame(
    {
        "Clean-trained": [
            np.mean(residuals_clean),
            np.std(residuals_clean),
            np.max(np.abs(residuals_clean)),
        ],
        "Noise-aware": [
            np.mean(residuals_aware),
            np.std(residuals_aware),
            np.max(np.abs(residuals_aware)),
        ],
    },
    index=[
        "Mean signed error (nm)",
        "Residual standard deviation (nm)",
        "Maximum absolute error (nm)",
    ],
)

print("\nResidual comparison")
print(residual_summary.to_string(float_format=lambda v: f"{v:.6f}"))


# Save results

results.to_csv("sans_krr_model_comparison.csv")
residual_summary.to_csv("sans_krr_residual_summary.csv")

print("\nFigures and result tables saved.")
```
