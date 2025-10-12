%% compare_iterative_MAP_MMSE_3D.m
% DEMO: Iterative MAP vs. MMSE volume alignment using an SO(3) template bank.
% Author: Amnon Balanov
% Last updated: 12.10.2025

%% ------------------------------ Configuration ----------------------------
cfg = struct();

%% Path handling (relative to this script)
paths.script_dir   = get_script_dir();
paths.project_root = fileparts(paths.script_dir);
paths.assets_dir   = fullfile(paths.project_root, 'assets');
paths.rules_dir    = fullfile(paths.project_root, 'sphere_rules');  % for quadrature rules

% Volumes
cfg.path.volX = fullfile(paths.assets_dir, 'emdb_2984.mat');     % adjust if needed
cfg.path.volY = fullfile(paths.assets_dir, 'S80_ribosome.mat');  % adjust if needed
cfg.field.volX = 'volume';        % expect volume(1,:,:,:)
cfg.field.volY = 'original_vol';  % expect original_vol(1,:,:,:)

% Preprocessing
cfg.d = 32;
cfg.norm_each = true;
cfg.pre_rotate.volX = struct('angle', 90,  'axis', [0 1 1]); % degrees
cfg.pre_rotate.volY = struct('angle', 300, 'axis', [0 0 1]); % degrees

% Observations
cfg.N = 3000;            % number of observations
cfg.sigma = 0.4;         % noise std. Start with lower noises.
cfg.max_iters = 100;
cfg.tol = 1e-6;

% Iterative template bank: choose RANDOM or QUADRATURE
cfg.use_quadrature = true;   

% If RANDOM grid:
cfg.L = 600;                 % number of templates if random

% If QUADRATURE grid:
cfg.quad.L         = 4;      % rule parameters (increase for denser grid)
cfg.quad.P         = 3;
cfg.quad.k         = 2;
cfg.quad.reduction = 0;

% Optional prior for template weights (applies to both random and quadrature)
cfg.template_prior.mode  = 'uniform';   % 'uniform' or 'iso-gaussian'
cfg.template_prior.sigma = 0.5;         % radians (if iso-gaussian)
cfg.template_prior.center.type = 'identity';
% e.g. axis-angle center:
% cfg.template_prior.center = struct('type','axis-angle','axis',[0 0 1],'angle',0.3);

% Rotations + templates
cfg.interp = 'cubic';
cfg.mode   = 'crop';

% RNG
cfg.seed = 1;

%% ------------------------------ Setup & Data -----------------------------
if ~isempty(cfg.seed), rng(cfg.seed); end

Sx = load(cfg.path.volX);
Sy = load(cfg.path.volY);

vol_x = Sx.(cfg.field.volX);  if ndims(vol_x)==4, vol_x = squeeze(vol_x(1,:,:,:)); end
vol_y = Sy.(cfg.field.volY);  if ndims(vol_y)==4, vol_y = squeeze(vol_y(1,:,:,:)); end

vol_x = imresize3(double(vol_x), [cfg.d cfg.d cfg.d]);
vol_y = imresize3(double(vol_y), [cfg.d cfg.d cfg.d]);

vol_x = imrotate3(vol_x, cfg.pre_rotate.volX.angle, cfg.pre_rotate.volX.axis, cfg.interp, cfg.mode);
vol_y = imrotate3(vol_y, cfg.pre_rotate.volY.angle, cfg.pre_rotate.volY.axis, cfg.interp, cfg.mode);

if cfg.norm_each
    vol_x = vol_x / max(norm(vol_x,'fro'), eps);
    vol_y = vol_y / max(norm(vol_y,'fro'), eps);
end

%% ---------------------- Ground-truth rotations & signals -----------------
% Observations remain RANDOM/Haar
[gt_axis, gt_angle] = generateSO3randomRotations(cfg.N);

sig = zeros(cfg.d, cfg.d, cfg.d, cfg.N, 'like', vol_y);
for i = 1:cfg.N
    y_clean = imrotate3(vol_y, rad2deg(gt_angle(i)), gt_axis(i,:), cfg.interp, cfg.mode);
    sig(:,:,:,i) = y_clean + cfg.sigma * randn(cfg.d, cfg.d, cfg.d, 'like', vol_y);
end

%% --------------------------- SO(3) template bank -------------------------
if cfg.use_quadrature
    % Deterministic SO(3) grid + quadrature weights
    quad_cfg = struct('L', cfg.quad.L, 'P', cfg.quad.P, 'k', cfg.quad.k, ...
                      'reduction', cfg.quad.reduction, 'rules_dir', paths.rules_dir);

    [grid_axis, grid_angle, grid_R, w] = generateSO3quadratureRotations(quad_cfg, cfg.template_prior);
else
    % Random (Haar) grid with optional prior reweighting of weights
    [grid_axis, grid_angle, grid_R, w] = generateSO3randomRotations(cfg.L, cfg.template_prior);
end

Lgrid = numel(w);
grid_Ri = zeros(3,3,Lgrid);
for l = 1:Lgrid
    grid_Ri(:,:,l) = axang2rotm([grid_axis(l,:), -grid_angle(l)]);
end

%% ---------------------------- Iterative Estimation -----------------------
x_MAP  = vol_x;     % hard-assignment state
x_MMSE = vol_x;     % soft-assignment state

for tt = 1:cfg.max_iters
    fprintf('Iteration %d/%d\n', tt, cfg.max_iters);

    % Build templates from current estimates at the chosen grid nodes
    volT_MAP  = zeros(cfg.d,cfg.d,cfg.d,Lgrid, 'like', vol_x);
    volT_MMSE = zeros(cfg.d,cfg.d,cfg.d,Lgrid, 'like', vol_x);
    parfor l = 1:Lgrid
        volT_MAP(:,:,:,l)  = imrotate3(x_MAP,  rad2deg(grid_angle(l)), grid_axis(l,:), cfg.interp, cfg.mode);
        volT_MMSE(:,:,:,l) = imrotate3(x_MMSE, rad2deg(grid_angle(l)), grid_axis(l,:), cfg.interp, cfg.mode);
    end

    x_MAP_next  = zeros(cfg.d,cfg.d,cfg.d, 'like', vol_x);
    x_MMSE_next = zeros(cfg.d,cfg.d,cfg.d, 'like', vol_x);

    parfor i = 1:cfg.N
        y = sig(:,:,:,i);

        % ------- MAP: pick the node with max inner product -------
        corr_map = squeeze(sum(y .* volT_MAP, [1 2 3]));     % Lgrid x 1
        [~, idx] = max(corr_map);
        y_map_aligned = imrotate3(y, -rad2deg(grid_angle(idx)), grid_axis(idx,:), cfg.interp, cfg.mode);
        x_MAP_next = x_MAP_next + y_map_aligned;

        % ------- MMSE: softmax over nodes WITH weights w -------
        corr_mmse = squeeze(sum(y .* volT_MMSE, [1 2 3]));   % Lgrid x 1
        exps = corr_mmse / (cfg.sigma^2);
        mexp = max(exps);
        num  = exp(exps - mexp);
        Z    = sum(w .* num);                                % partition with node weights
        p    = (w .* num) / max(Z, eps);                     % posterior over nodes (sum=1)

        % Weighted rotation average → project to SO(3)
        M = sum( grid_R .* reshape(p,1,1,[]), 3 );
        R_mmse = project_to_SO3(M);
        ax = rotm2axang(R_mmse);

        y_mmse_aligned = imrotate3(y, -rad2deg(ax(4)), ax(1:3), cfg.interp, cfg.mode);
        x_MMSE_next = x_MMSE_next + y_mmse_aligned;
    end

    % Averages
    x_MAP_new  = x_MAP_next  / cfg.N;
    x_MMSE_new = x_MMSE_next / cfg.N;

    % Convergence check
    dMAP  = norm(x_MAP - x_MAP_new,   'fro');
    dMMSE = norm(x_MMSE - x_MMSE_new, 'fro');
    fprintf('  ΔMAP = %.3g,  ΔMMSE = %.3g\n', dMAP, dMMSE);

    x_MAP  = x_MAP_new;
    x_MMSE = x_MMSE_new;

    if dMAP < cfg.tol && dMMSE < cfg.tol
        fprintf('Converged at iteration %d.\n', tt);
        break;
    end
end

%% ----------------------------- Final Alignment ---------------------------
[~, x_MAP_aligned]  = findBest3Dalignment(vol_y, x_MAP,  struct('L', 2000, 'interp', cfg.interp, 'mode', cfg.mode));
[~, x_MMSE_aligned] = findBest3Dalignment(vol_y, x_MMSE, struct('L', 2000, 'interp', cfg.interp, 'mode', cfg.mode));

%% ------------------------------- Visualization ---------------------------
% Show 3D volumes in a 2x2 dashboard with clear titles.
% Panels:
%   (1) MAP-aligned estimate          (x_MAP_aligned)
%   (2) MMSE-aligned estimate         (x_MMSE_aligned)
%   (3) Initial / reference volume X  (vol_x)
%   (4) Target volume Y               (vol_y)

% 1) Prepare display volumes: robust scaling to uint8
V1 = scale_volume_to_uint8(x_MAP_aligned);   % MAP-aligned reconstruction
V2 = scale_volume_to_uint8(x_MMSE_aligned);  % MMSE-aligned reconstruction
V3 = scale_volume_to_uint8(vol_x);           % initial/reference X
V4 = scale_volume_to_uint8(vol_y);           % target Y

% 2) UI figure + grid with titled panels
fig = uifigure('Name','Iterative MAP vs. MMSE — Volume Viewer');

gl = uigridlayout(fig,[2,2]);
gl.RowHeight     = {'1x','1x'};
gl.ColumnWidth   = {'1x','1x'};
gl.Padding       = [8 8 8 8];
gl.RowSpacing    = 6;
gl.ColumnSpacing = 6;

p1 = uipanel(gl, 'Title','MAP-aligned estimate (x\_MAP\_aligned)');
p2 = uipanel(gl, 'Title','MMSE-aligned estimate (x\_MMSE\_aligned)');
p3 = uipanel(gl, 'Title','Initial/reference X (vol\_x)');
p4 = uipanel(gl, 'Title','Target Y (vol\_y)');

% 3) Create viewer3d in each panel (set background on the viewer, not on volshow)
bg = [0 0 0];                         % black background
vwr1 = viewer3d('Parent', p1, 'BackgroundColor', bg);
vwr2 = viewer3d('Parent', p2, 'BackgroundColor', bg);
vwr3 = viewer3d('Parent', p3, 'BackgroundColor', bg);
vwr4 = viewer3d('Parent', p4, 'BackgroundColor', bg);

% 4) Consistent volume styling
cm  = gray(256);                      % colormap
alp = linspace(0,1,256).^0.8;         % slightly nonlinear opacity

% Use 'RenderingStyle' (new API), DO NOT pass BackgroundColor to volshow
h1 = volshow(V1*0.005, 'Parent', vwr1, 'Colormap', cm, 'Alphamap', alp, 'RenderingStyle','VolumeRendering');
h2 = volshow(V2*0.005, 'Parent', vwr2, 'Colormap', cm, 'Alphamap', alp, 'RenderingStyle','VolumeRendering');
h3 = volshow(V3*0.005, 'Parent', vwr3, 'Colormap', cm, 'Alphamap', alp, 'RenderingStyle','VolumeRendering');
h4 = volshow(V4*0.005, 'Parent', vwr4, 'Colormap', cm, 'Alphamap', alp, 'RenderingStyle','VolumeRendering');

% 5) (Optional) synchronize camera across viewers for easier comparison
try
    vwr1.CameraPosition = vwr4.CameraPosition;
    vwr1.CameraTarget   = vwr4.CameraTarget;
    vwr1.CameraUpVector = vwr4.CameraUpVector;

    vwr2.CameraPosition = vwr4.CameraPosition;
    vwr2.CameraTarget   = vwr4.CameraTarget;
    vwr2.CameraUpVector = vwr4.CameraUpVector;

    vwr3.CameraPosition = vwr4.CameraPosition;
    vwr3.CameraTarget   = vwr4.CameraTarget;
    vwr3.CameraUpVector = vwr4.CameraUpVector;
catch
    
end


%% ------------------------------- Utilities -------------------------------
function R = project_to_SO3(M)
%PROJECT_TO_SO3 Closest rotation via SVD (orthogonal Procrustes with det=+1).
    [U,~,V] = svd(M);
    S = eye(3);
    if det(U*V') < 0, S(3,3) = -1; end
    R = U*S*V';
end

function d = get_script_dir()
%GET_SCRIPT_DIR Directory of the current script, even if run as a function.
    fp = mfilename('fullpath');
    if isempty(fp), d = pwd; else, d = fileparts(fp); end
end

function U8 = scale_volume_to_uint8(V)
%SCALE_VOLUME_TO_UINT8 Robustly map a 3D volume to uint8 for volshow.
% Clips to the 2nd–98th percentile to avoid outliers flattening contrast.

    V  = double(V);
    lo = prctile(V(:), 2);
    hi = prctile(V(:), 98);
    if ~isfinite(lo) || ~isfinite(hi) || hi <= lo
        lo = min(V(:)); hi = max(V(:));
        if hi <= lo, U8 = uint8(zeros(size(V))); return; end
    end
    Vc = min(max(V, lo), hi);
    U8 = uint8(255 * (Vc - lo) / (hi - lo));
end
