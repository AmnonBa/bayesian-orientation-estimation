bayesian-orientation-estimation

Code accompanying the paper “A Bayesian Perspective on Orientation Determination in Cryo-EM with Applications to Structural Heterogeneity Analysis.”

══════════════════════════════════════════════════════════════════════════════
SO(3) ORIENTATION ESTIMATION & ITERATIVE ALIGNMENT — MATLAB TOOLKIT
══════════════════════════════════════════════════════════════════════════════

Demos and utilities for MAP/MMSE orientation estimation on SO(3),
supporting both RANDOM (Haar) and QUADRATURE grids. Includes iterative
alignment and clean 3D visualization presets.

──────────────────────────────────────────────────────────────────────────────
▌ CONTENTS
──────────────────────────────────────────────────────────────────────────────
main/
├─ compare_iterative_MAP_MMSE_3D.m ← Iterative MAP vs. MMSE alignment (EM vs hard-assignment)
├─ compare_so3_MAP_vs_MMSE.m ← MAP/MMSE error vs. SNR sweep
├─ generateSO3randomRotations.m ← Random (Haar / iso-Gaussian) SO(3) nodes
├─ generateSO3quadratureRotations.m ← Quadrature SO(3) nodes + weights
└─ findBest3Dalignment.m ← Random alignment search (NCC / dot)

assets/
├─ emdb_2984.mat (field: volume)
└─ S80_ribosome.mat (field: original_vol)

quadrature/ ← all files realted to the generation of the SO(3) quadrature points.
├─sphere_rules/ ←  *.dat rule files
results/ ← auto-created for plots & figures

──────────────────────────────────────────────────────────────────────────────
▌ REQUIREMENTS
──────────────────────────────────────────────────────────────────────────────
• MATLAB R2021a+ (tested)
• Image Processing Toolbox (imrotate3, volshow, viewer3d)
• A loader:
get_SO3_rule(L, P, k, reduction, ...) → [alpha beta gamma w] (radians)
Place it in quadrature/, and any rule files (*.dat) in sphere_rules/.

Visual API notes:
• Use viewer3d for background color; volshow uses 'RenderingStyle' (no old “Renderer”).
• Do NOT set BackgroundColor on volshow’s Volume object.

──────────────────────────────────────────────────────────────────────────────
▌ QUICK START
──────────────────────────────────────────────────────────────────────────────

Add paths & enter project:
─────────────────────────────────────────────────────────
addpath('main','quadrature','sphere_rules'); % as needed
cd main
─────────────────────────────────────────────────────────

Run the iterative alignment demo:
─────────────────────────────────────────────────────────
compare_iterative_MAP_MMSE_3D
─────────────────────────────────────────────────────────
• Loads assets/emdb_2984.mat (X) & assets/S80_ribosome.mat (Y)
• Generates N noisy rotated observations of Y
• Builds SO(3) template bank (RANDOM or QUADRATURE)
• Iteratively refines MAP (hard) and MMSE (soft) reconstructions
• Shows a 2×2 dashboard: MAP, MMSE, reference X, target Y

SNR sweep (MAP vs. MMSE orientation error):
─────────────────────────────────────────────────────────
compare_so3_MAP_vs_MMSE
─────────────────────────────────────────────────────────
• Sweeps cfg.snr_list, computes mean geodesic error (radians)
• Saves plot to results/ if configured

──────────────────────────────────────────────────────────────────────────────
▌ KEY PARAMETERS (CHEAT SHEET)
──────────────────────────────────────────────────────────────────────────────
General (both demos)
• cfg.d : volume side length after resize (e.g., 32)
• cfg.interp : 'cubic' (recommended)
• cfg.mode : 'crop' (keep size)

Iterative demo (compare_iterative_MAP_MMSE_3D.m)
• cfg.N : number of observations (e.g., 3000)
• cfg.sigma : noise std (start modest; e.g., 4e-5 or 0.4 for stress)
• cfg.max_iters : outer iterations
• cfg.tol : stopping tolerance
• cfg.use_quadrature : true → QUAD nodes; false → RANDOM nodes

Random grid:
• cfg.L : number of nodes/templates

Quadrature grid:
• cfg.quad.L, cfg.quad.P, cfg.quad.k, cfg.quad.reduction
↑ Increase for denser coverage (cost ↑)

Optional template PRIOR (applies to weights w, not node placement):
• cfg.template_prior.mode : 'uniform' | 'iso-gaussian'
• cfg.template_prior.sigma : std in radians (iso-Gaussian)
• cfg.template_prior.center.type : 'identity' | 'axis-angle' | 'matrix'
(if 'axis-angle', set .axis and .angle; if 'matrix', set .R0)

SNR sweep (compare_so3_MAP_vs_MMSE.m)
• cfg.N_obs : observations per SNR
• cfg.snr_list : vector, e.g., logspace(-5, -1, 12)
• cfg.use_quadrature : as above
• cfg.random_grid_L or cfg.quad.* : as above

──────────────────────────────────────────────────────────────────────────────
▌ HOW THINGS WORK
──────────────────────────────────────────────────────────────────────────────
MAP (hard assignment):

Build template bank {Rℓ} on SO(3) (RANDOM or QUADRATURE).

For each observation yᵢ, choose ℓ* = argmax_ℓ <yᵢ, Tℓ>.

Back-rotate yᵢ by Rℓ*⁻¹ and average over i.

MMSE (soft assignment):

Compute correlations cℓ = <yᵢ, Tℓ> and softmax exps(cℓ/σ²).

Include node weights wℓ (Haar or prior-reweighted): pℓ x wℓ·softmax(cℓ/σ²).

Average rotations M = Σℓ pℓ Rℓ, then project M → SO(3) via SVD (det=+1).

Back-rotate yᵢ by the MMSE rotation and average over i.

Quadrature prior reweighting (iso-Gaussian about R₀):
• wℓ ← wℓ · exp(−θ(Rℓ,R₀)² / (2σ²)), then renormalize Σℓ wℓ = 1.

──────────────────────────────────────────────────────────────────────────────
▌ APIs YOU’LL CALL
──────────────────────────────────────────────────────────────────────────────
Random nodes:
[axis, angle, R, w] = generateSO3randomRotations(N, prior_cfg)
• prior_cfg.mode = 'uniform' → Haar via unit quaternions
• prior_cfg.mode = 'iso-gaussian' → axis, S²; angle from geodesic Gaussian
• w returns equal weights (or reweighted if you extend it)

Quadrature nodes:
[axis, angle, R, w] = generateSO3quadratureRotations(quad_cfg, prior_cfg)
• Uses get_SO3_rule to produce ZYZ angles + Haar weights
• Optional iso-Gaussian reweighting of w

Random alignment search:
[axang, vol2_best, idx, score] = findBest3Dalignment(vol1, vol2, opts)
• opts.normalize = true → NCC; falls back to dot product if degenerate

──────────────────────────────────────────────────────────────────────────────
▌ VISUALIZATION (NEW VIEWER API)
──────────────────────────────────────────────────────────────────────────────
• Background color: set on viewer3d, not on volshow
• Use RenderingStyle with volshow

Example (from the demo’s 2×2 dashboard):
vwr = viewer3d('Parent', panel, 'BackgroundColor',[0 0 0]);
volshow(data_uint8, 'Parent', vwr, ...
'Colormap', gray(256), ...
'Alphamap', linspace(0,1,256).^0.8, ...
'RenderingStyle','VolumeRendering');

Volumes are auto-scaled to uint8 with robust percentile clipping (2–98%).

──────────────────────────────────────────────────────────────────────────────
▌ TROUBLESHOOTING
──────────────────────────────────────────────────────────────────────────────
• “BackgroundColor is not supported on the Volume object.”
→ Set the background on viewer3d, not on volshow.

• “Renderer property is no longer supported.”
→ Use 'RenderingStyle' on volshow; do not set obsolete 'Renderer'.

• get_SO3_rule missing or no rule files
→ Ensure quadrature/ is on MATLAB path; put rule files in sphere_rules/.

• Blank/low-contrast volumes
→ Check normalization; verify robust scaling; try a larger cfg.d and
a gentler opacity map (e.g., linspace(0,1,256).^0.6).

• Slow/huge memory usage
→ Reduce cfg.N, cfg.L (or quadrature density), cfg.max_iters; start with
higher σ to test the pipeline before scaling up.

──────────────────────────────────────────────────────────────────────────────
▌ CITATION
──────────────────────────────────────────────────────────────────────────────
If you use this code, please cite:
title={Bayesian Perspective for Orientation Estimation in Cryo-EM and Cryo-ET},
author={Xu, Sheng and Balanov, Amnon and Singer, Amit and Bendory, Tamir},
journal={bioRxiv},
pages={2025--10},
year={2025},
publisher={Cold Spring Harbor Laboratory}

