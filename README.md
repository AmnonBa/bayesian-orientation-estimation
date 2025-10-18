# Bayesian Orientation Estimation on SO(3)

This MATLAB toolkit accompanies the paper:  
**“A Bayesian Perspective on Orientation Determination in Cryo-EM with Applications to Structural Heterogeneity Analysis.”**

It provides implementations and demos for **MAP** (Maximum A Posteriori) and **MMSE** (Minimum Mean-Square Error) orientation estimation on **SO(3)**, including iterative alignment and visualization tools.  
Both **random (Haar or isotropic Gaussian)** and **quadrature-based** SO(3) grids are supported.

---

## 🧭 Overview

This package enables Bayesian orientation estimation for 3D Cryo-EM and Cryo-ET volumes. It includes:

- MAP and MMSE estimation on SO(3)
- Random and quadrature SO(3) grids
- Iterative alignment and EM-like updates
- Clean 3D visualization presets for MATLAB (`viewer3d`, `volshow`)
- Automated benchmarking for orientation error vs. SNR

---

## 📂 Repository Structure

```
bayesian-orientation-estimation/
│
├── main/
│   ├── compare_iterative_MAP_MMSE_3D.m      # Iterative MAP vs. MMSE alignment
│   ├── compare_so3_MAP_vs_MMSE.m            # SNR sweep (MAP vs. MMSE error)
│   ├── generateSO3randomRotations.m         # Random SO(3) nodes (Haar / iso-Gaussian)
│   ├── generateSO3quadratureRotations.m     # Quadrature SO(3) nodes + weights
│   └── findBest3Dalignment.m                # Random 3D alignment (NCC / dot product)
│
├── quadrature/                              # Generation of SO(3) quadrature points
│   └── sphere_rules/                        # *.dat rule files for quadrature
│
├── assets/                                  # Example 3D volumes
│   ├── emdb_2984.mat                        # Reference volume
│   └── S80_ribosome.mat                     # Target volume
│
├── results/                                 # Auto-created for plots and figures
└── README.md
```

---

## ⚙️ Requirements

- **MATLAB R2021a** or later (tested)
- **Image Processing Toolbox** (for `imrotate3`, `volshow`, `viewer3d`)
- SO(3) quadrature loader function:  
  ```matlab
  [alpha, beta, gamma, w] = get_SO3_rule(L, P, k, reduction, ...);
  ```
  Place it in `quadrature/`, with rule files (`*.dat`) in `sphere_rules/`.

**Visualization notes:**
- Use `viewer3d` for background color customization.
- `volshow` uses the `'RenderingStyle'` parameter (no `Renderer` property).
- Avoid setting `BackgroundColor` on the `Volume` object directly.

---

## 🚀 Quick Start

### 1️⃣ Setup
Add required paths and move into the project directory:
```matlab
addpath('main','quadrature','quadrature/sphere_rules');
cd main
```

### 2️⃣ Iterative Alignment Demo
```matlab
compare_iterative_MAP_MMSE_3D
```
This demo:
- Loads `assets/emdb_2984.mat` (reference X) and `assets/S80_ribosome.mat` (target Y)
- Generates N noisy rotated observations of Y
- Builds an SO(3) template bank (RANDOM or QUADRATURE)
- Performs iterative MAP (hard) and MMSE (soft) reconstructions
- Displays a 2×2 dashboard: MAP, MMSE, reference X, and target Y

### 3️⃣ SNR Sweep (MAP vs. MMSE Orientation Error)
```matlab
compare_so3_MAP_vs_MMSE
```
This demo:
- Sweeps over `cfg.snr_list`
- Computes mean geodesic error (in radians)
- Optionally saves plots to the `results/` folder

---

## 🧠 Key Parameters

### General
| Parameter | Description | Example |
|------------|-------------|----------|
| `cfg.d` | Volume side length after resize | 32 |
| `cfg.interp` | Interpolation mode | `'cubic'` |
| `cfg.mode` | Resizing behavior | `'crop'` |

### Iterative Demo (`compare_iterative_MAP_MMSE_3D.m`)
| Parameter | Description | Example |
|------------|-------------|----------|
| `cfg.N` | Number of observations | 3000 |
| `cfg.sigma` | Noise standard deviation | 0.1 |
| `cfg.max_iters` | Number of outer iterations | 10 |
| `cfg.tol` | Stopping tolerance | 1e-4 |
| `cfg.use_quadrature` | Use quadrature grid instead of random | `true` |

**Random grid:**  
`cfg.L` — number of random nodes/templates

**Quadrature grid:**  
`cfg.quad.L`, `cfg.quad.P`, `cfg.quad.k`, `cfg.quad.reduction` — controls density and accuracy

### Optional Template Prior
| Parameter | Description |
|------------|-------------|
| `cfg.template_prior.mode` | `'uniform'` or `'iso-gaussian'` |
| `cfg.template_prior.sigma` | Standard deviation (radians) for Gaussian prior |
| `cfg.template_prior.center` | `'identity'`, `'axis-angle'`, or `'matrix'` |

### SNR Sweep (`compare_so3_MAP_vs_MMSE.m`)
| Parameter | Description |
|------------|-------------|
| `cfg.N_obs` | Observations per SNR |
| `cfg.snr_list` | Vector of SNR values (e.g. `logspace(-5,-1,12)`) |
| `cfg.use_quadrature` | Whether to use quadrature nodes |
| `cfg.random_grid_L` / `cfg.quad.*` | Node definitions |

---

## ⚙️ Algorithmic Summary

### MAP (Hard Assignment)
1. Build SO(3) template bank `{R_ℓ}` (RANDOM or QUADRATURE).
2. For each observation `yᵢ`, choose:
   ```math
   ℓ* = argmax_ℓ ⟨yᵢ, T_ℓ⟩
   ```
3. Back-rotate `yᵢ` by `R_ℓ*^{-1}` and average over all `i`.

### MMSE (Soft Assignment)
1. Compute correlations `c_ℓ = ⟨yᵢ, T_ℓ⟩`.
2. Apply softmax: `p_ℓ ∝ w_ℓ exp(c_ℓ / σ²)`.
3. Average rotations `M = Σℓ p_ℓ R_ℓ` and project `M → SO(3)` via SVD (det=+1).
4. Back-rotate each `yᵢ` by the MMSE rotation and average.

### Quadrature Prior Reweighting
For an isotropic Gaussian prior centered at `R₀`:
```math
w_ℓ ← w_ℓ · exp(-θ(R_ℓ, R₀)² / (2σ²)), \quad Σ_ℓ w_ℓ = 1
```

---

## 📘 API Reference

### Random Nodes
```matlab
[axis, angle, R, w] = generateSO3randomRotations(N, prior_cfg)
```
- `prior_cfg.mode = 'uniform'` → Haar via unit quaternions  
- `prior_cfg.mode = 'iso-gaussian'` → Axis sampled from S², angle from Gaussian  
- Returns equal or reweighted node weights `w`

### Quadrature Nodes
```matlab
[axis, angle, R, w] = generateSO3quadratureRotations(quad_cfg, prior_cfg)
```
- Uses `get_SO3_rule()` to load quadrature rules (`ZYZ` angles + weights)  
- Supports optional Gaussian reweighting of `w`

### Random Alignment
```matlab
[axang, vol2_best, idx, score] = findBest3Dalignment(vol1, vol2, opts)
```
- `opts.normalize = true` → normalized cross-correlation (NCC)  
- Falls back to dot product if normalization fails

---

## 📚 Citation

If you use this code, please cite:

> Xu, Sheng; Balanov, Amnon; Singer, Amit; Bendory, Tamir.  
> *"A Bayesian Perspective for Orientation Estimation in Cryo-EM and Cryo-ET."*  
> **bioRxiv**, 2025, Cold Spring Harbor Laboratory.

---

## 👤 Authors

**Sheng Xu**, **Amnon Balanov**, **Amit Singer**, and **Tamir Bendory**  
Department of Electrical Engineering, Tel Aviv University

---

## 🧾 License

This software is provided for academic and research use only.  
© 2025 The Authors. All rights reserved.
