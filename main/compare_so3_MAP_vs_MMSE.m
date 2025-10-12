%% compare_so3_MAP_vs_MMSE.m
% Compare MAP vs. MMSE orientation estimation on SO(3) using a template bank.
% Author: Amnon Balanov
% Last updated: 12.10.2025
%
% Overview:
%   1) Load and normalize a 3D volume
%   2) Build an SO(3) template bank (Monte-Carlo or quadrature-based)
%   3) Generate ground-truth rotations and noisy observations for multiple SNRs
%   4) Estimate rotations via MAP and MMSE (softmax-weighted average projected to SO(3))
%   5) Report geodesic errors and plot mean error vs. SNR
%
% Requirements:
%   - MATLAB R2021a+ (tested), Image Processing Toolbox (for imrotate3)
%   - If using quadrature: sphere_rules/*.dat as in your original setup
%
% Notes:
%   - Templates are re-normalized to Frobenius norm 1 to match the likelihood model.
%   - The MMSE density uses weights from the SO(3) rule when enabled.

%% --------------------------- Path Handling -------------------------------
% Resolve everything relative to this script’s folder, so it works no matter
% where you run it from.
paths.script_dir     = get_script_dir();
paths.project_root   = fileparts(paths.script_dir);
paths.assets_dir     = fullfile(paths.project_root, 'assets');
paths.rules_dir      = fullfile(paths.project_root, 'sphere_rules');
paths.results_dir    = fullfile(paths.project_root, 'results');

% Make sure output dirs exist (assets/rules are read-only; results is write)
ensure_dir(paths.results_dir);

%% --------------------------- Configuration ------------------------------
cfg = struct();

% Data / volume
% You can set an absolute path, or a relative path (relative to the script).
% If empty, default to assets/S80_ribosome.mat.
cfg.volume_mat_path   = fullfile('..','assets','S80_ribosome.mat');  % relative example
cfg.volume_field_name = 'original_vol';   % variable inside the MAT file
cfg.take_index_first  = true;             % original_vol(1,:,:,:) slice if true
cfg.d                 = 32;               % side length after resize

% Template bank (SO(3) grid)
cfg.use_quadrature    = false;             % if false: use random grid

% SO(3) quadrature knobs
cfg.quad.L            = 6;
cfg.quad.P            = 5;
cfg.quad.k            = 4;
cfg.quad.reduction    = 0;
cfg.random_grid_L     = 10250;              % number of templates if not using quadrature

% Template prior for the TEMPLATE BANK (applies to quadrature & random)
% 'uniform'           : keep weights as-is (Haar for quadrature; equal for random)
% 'iso-gaussian'      : reweight by exp( -theta(R,R0)^2 / (2*sigma^2) )
cfg.template_prior.mode  = 'uniform';      % 'uniform' or 'iso-gaussian'
cfg.template_prior.sigma = 0.1;            % radians; std of geodesic angle for 'iso-gaussian'
% Center of the prior (mean rotation R0):
cfg.template_prior.center.type = 'identity';  % 'identity' | 'axis-angle' | 'matrix'

% If using 'axis-angle', also set:
% cfg.template_prior.center.axis  = [0 0 1];
% cfg.template_prior.center.angle = 0.3;    % radians
% If using 'matrix', also set:
% cfg.template_prior.center.R0    = eye(3);

% Experiment design
cfg.N_obs             = 2000;             % number of observations per point
cfg.snr_list          = logspace(-5, -1, 12); % SNR per-sample (signal power / noise power)
cfg.rng_seed          = 42;               % reproducibility

% Rotation prior for ground truth
cfg.gt.mode           = 'iso-gaussian';   % 'uniform' or 'iso-gaussian'
cfg.gt.sigma          = 0.1;              % only used if 'iso-gaussian'

% Implementation details
cfg.rotate_interp     = 'cubic';          % imrotate3 interpolation
cfg.rotate_mode       = 'crop';           % keep size dxdxd
cfg.normalize_templates = true;           % enforce ||template||_F = 1

% Plotting
cfg.show_plot         = true;
% If empty: no save.
% If a directory: will save as <dir>/results_mean_geo_error.png
% If a filename (with extension): will save exactly there (dir auto-created).
cfg.save_plot_path    = fullfile(paths.results_dir); % directory example

% Quadrature rules location (folder with your *.dat rule files)
cfg.rules_dir         = paths.rules_dir;

%% ------------------------- Reproducibility ------------------------------
rng(cfg.rng_seed);

%% ----------------------------- Resolve Paths -----------------------------
% Volume path: absolute if already absolute; otherwise relative to script
if isempty(cfg.volume_mat_path)
    candidate = fullfile(paths.assets_dir, 'S80_ribosome.mat');
else
    candidate = resolve_path(cfg.volume_mat_path, paths.script_dir);
end
assert(exist(candidate,'file')==2, ...
    'Volume MAT not found. Tried: %s', candidate);
cfg.volume_mat_path = candidate;

% Save path: allow "", directory, or filename
if ~isempty(cfg.save_plot_path)
    cfg.save_plot_path = resolve_path(cfg.save_plot_path, paths.script_dir);
    if isfolder(cfg.save_plot_path)
        % given a directory -> pick a default filename there
        cfg.save_plot_path = fullfile(cfg.save_plot_path, 'results_mean_geo_error.png');
    else
        % given a file path -> ensure its folder exists
        ensure_dir(fileparts(cfg.save_plot_path));
    end
end


%% ----------------------------- Load volume ------------------------------
S = load(cfg.volume_mat_path);
assert(isfield(S, cfg.volume_field_name), ...
    'Field "%s" not found in MAT file: %s', cfg.volume_field_name, cfg.volume_mat_path);

vol_full = S.(cfg.volume_field_name);

% Pick a volume: many files store a stack; we take the first if requested
if cfg.take_index_first
    % Expect vol_full(1,:,:,:) with size [1 d0 d0 d0]
    vol = squeeze(vol_full(1,:,:,:));
else
    vol = vol_full;
end

% Resize and normalize
d = cfg.d;
vol = imresize3(double(vol), [d d d]);

% Normalize to Frobenius norm 1 (signal power = 1)
vol = vol / max(norm(vol, 'fro'), eps);

%% -------- SO(3) template bank: rotations, weights, and templates ---------
if cfg.use_quadrature
    quad_cfg = struct('L', cfg.quad.L, 'P', cfg.quad.P, 'k', cfg.quad.k, 'reduction', cfg.quad.reduction);
    % If your loader needs it: quad_cfg.rules_dir = cfg.rules_dir;
    [grid_axis, grid_angle, grid_R, w] = generateSO3quadratureRotations(quad_cfg, cfg.template_prior);
else
    [grid_axis, grid_angle, grid_R, w] = generateSO3randomRotations(cfg.random_grid_L, cfg.template_prior);
end
Lgrid = numel(w);

% Precompute rotated templates
vol_templates = zeros(d,d,d,Lgrid, 'like', vol);
for i = 1:Lgrid
    vol_i = imrotate3(vol, rad2deg(grid_angle(i)), grid_axis(i,:), cfg.rotate_interp, cfg.rotate_mode);
    if cfg.normalize_templates
        nrm = norm(vol_i, 'fro');
        if nrm > 0
            vol_i = vol_i / nrm;
        end
    end
    vol_templates(:,:,:,i) = vol_i;
end


%% -------------- Ground truth rotations and noisy observations -----------
N = cfg.N_obs;
[gt_axis, gt_angle] = generateSO3randomRotations(N, cfg.gt);
gt_R = zeros(N,3,3);
for i = 1:N
    gt_R(i,:,:) = axang2rotm([gt_axis(i,:), gt_angle(i)]);
end

% Signal "power" under our normalization is ~1 (Fro norm 1).
% For each SNR value, sigma = sqrt(signal_var / SNR).
snr_list = cfg.snr_list;
sigma_list = sqrt( var(vol(:)) ./ snr_list );

%% ------------------ Storage for per-SNR evaluation -----------------------
R_MAP_geo    = zeros(N, numel(snr_list));
R_MMSE_geo   = zeros(N, numel(snr_list));
R_MAP_mse    = zeros(N, numel(snr_list));
R_MMSE_mse   = zeros(N, numel(snr_list));

%% --------------------------- Main evaluation -----------------------------
for t = 1:numel(snr_list)
    fprintf('SNR %g (%d/%d)\n', snr_list(t), t, numel(snr_list));
    sigma = sigma_list(t);

    for i = 1:N
        % Generate noisy observation y_i
        y_clean = imrotate3(vol, rad2deg(gt_angle(i)), gt_axis(i,:), cfg.rotate_interp, cfg.rotate_mode);
        y = y_clean + sigma * randn(d, d, d);

        % Compute inner products <y, T_l>
        ip = sum(y .* vol_templates, [1 2 3]);           % 1x1x1xLgrid
        ip = reshape(ip, [Lgrid, 1]);                    % L x 1

        %% MAP estimate
        [~, idx_max] = max(ip.*w);
        R_map = grid_R(:,:,idx_max);

        R_MAP_mse(i,t)  = norm(squeeze(R_map) - squeeze(gt_R(i,:,:)), 'fro');
        R_MAP_geo(i,t)  = so3_geodesic_angle(R_map, squeeze(gt_R(i,:,:)));

        %% MMSE estimate (softmax over templates with weights)
        exponents = ip / (sigma^2);

        % Numerically stable softmax with quadrature weights
        mexp = max(exponents);
        num  = exp(exponents - mexp);            % L x 1
        Z    = sum(w .* num);                    % scalar
        dens = num / max(Z, eps);                % normalized densities

        % Posterior weights over templates:
        p = w .* dens;                           % L x 1, sums to 1

        % Weighted average of rotation matrices (in R^{3x3}); then project to SO(3)
        M = sum( grid_R .* reshape(p, 1,1,[]), 3 );    % 3x3
        R_mmse = project_to_SO3(M);

        R_MMSE_mse(i,t) = norm(R_mmse - squeeze(gt_R(i,:,:)), 'fro');
        R_MMSE_geo(i,t) = so3_geodesic_angle(R_mmse, squeeze(gt_R(i,:,:)));
    end
end

%% --------------------------------- Plot ----------------------------------
if cfg.show_plot
    figure; hold on;
    plot(snr_list, mean(R_MAP_geo, 1), '-o', 'DisplayName', 'MAP (mean geodesic)');
    plot(snr_list, mean(R_MMSE_geo,1), '-x', 'DisplayName', 'MMSE (mean geodesic)');
    set(gca, 'XScale', 'log'); grid on;
    xlabel('SNR'); ylabel('Geodesic error (radians)');
    title('SO(3) orientation estimation: MAP vs. MMSE');
    legend('Location','northeast');

    if ~isempty(cfg.save_plot_path)
        exportgraphics(gca, cfg.save_plot_path, 'Resolution', 200);
        fprintf('Saved plot to: %s\n', cfg.save_plot_path);
    end
end

%% ------------------------------ Utilities --------------------------------
function R = project_to_SO3(M)
%PROJECT_TO_SO3 Projects a 3x3 matrix M to the closest rotation (Procrustes).
    [U,~,V] = svd(M);
    S = eye(3);
    if det(U*V') < 0
        S(3,3) = -1; % enforce det=+1
    end
    R = U*S*V';
end

function th = so3_geodesic_angle(R1, R2)
%SO3_GEODESIC_ANGLE Geodesic distance on SO(3) between R1,R2 (both 3x3).
    C = R1' * R2;
    x = (trace(C) - 1) / 2;
    x = max(min(x, 1), -1); % clamp for numerical safety
    th = acos(x);
end

function out = resolve_path(p, base_dir)
%RESOLVE_PATH Return absolute path. If p is absolute, return as-is;
% otherwise resolve relative to base_dir.
    if isempty(p)
        out = p;
        return;
    end
    if is_absolute_path(p)
        out = p;
    else
        out = fullfile(base_dir, p);
    end
    % Normalize path (convert to absolute canonical form)
    out = char(java.io.File(out).getCanonicalPath());
end

function tf = is_absolute_path(p)
%IS_ABSOLUTE_PATH True if p is absolute on Windows or UNIX.
    if ispc
        % Drive-letter or UNC paths
        tf = (~isempty(regexp(p, '^[A-Za-z]:[\\/]', 'once'))) || startsWith(p,'\\') || startsWith(p,'//');
    else
        tf = startsWith(p, '/');
    end
end

function ensure_dir(d)
%ENSURE_DIR Create directory if it doesn't exist (no error if it already does).
    if isempty(d); return; end
    if exist(d,'dir') ~= 7
        mkdir(d);
    end
end

function d = get_script_dir()
%GET_SCRIPT_DIR Directory of the current script, even if run as a function.
    fp = mfilename('fullpath');
    if isempty(fp)
        % Fallback when running as a live script; use pwd
        d = pwd;
    else
        d = fileparts(fp);
    end
end


