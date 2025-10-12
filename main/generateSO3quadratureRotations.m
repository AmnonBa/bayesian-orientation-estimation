function [grid_axis, grid_angle, grid_R, w] = generateSO3quadratureRotations(quadrature_cfg, prior_cfg)
% Build an SO(3) node set from a quadrature rule.
%
% [axis, angle, R, w] = generateSO3quadratureRotations(quadrature_cfg, prior_cfg)
%
% Inputs
%   quadrature_cfg : struct with fields
%                    .L, .P, .k, .reduction  (rule parameters)
%                    (Optionally you can add .rules_dir if your loader needs it)
%   prior_cfg      : struct with fields
%                    .mode   = 'uniform' or 'iso-gaussian'
%                    .sigma  = std (radians) for geodesic Gaussian (if 'iso-gaussian')
%                    .center = see generateSO3randomRotations()
%
% Outputs (same API as random)
%   grid_axis  : L x 3 unit vectors
%   grid_angle : L x 1 angles in [0, pi]
%   grid_R     : 3 x 3 x L rotation matrices
%   w          : L x 1 quadrature weights (normalized), possibly reweighted by prior
%
% Notes
% - Assumes availability of get_SO3_rule(L,P,k,reduction, ...) that returns [alpha beta gamma w].
% - We reweight Haar quadrature weights by the prior density and renormalize.

L = quadrature_cfg.L;
P = quadrature_cfg.P;
k = quadrature_cfg.k;
r = quadrature_cfg.reduction;

% If your get_SO3_rule supports rules_dir, pass it here:
% SO3_rule = get_SO3_rule(L,P,k,r, quadrature_cfg.rules_dir);
SO3_rule = get_SO3_rule(L,P,k,r);

alpha = SO3_rule(:,1);
beta  = SO3_rule(:,2);
gamma = SO3_rule(:,3);
w     = SO3_rule(:,4);

Lgrid = numel(w);
grid_R     = zeros(3,3,Lgrid);
grid_axis  = zeros(Lgrid,3);
grid_angle = zeros(Lgrid,1);

for i = 1:Lgrid
    R = eul2rotm([alpha(i) beta(i) gamma(i)], 'ZYZ'); % radians
    grid_R(:,:,i) = R;
    axang = rotm2axang(R);
    grid_axis(i,:)  = axang(1:3);
    grid_angle(i,1) = axang(4);
end

% Ensure nonnegative base weights and normalize (Haar)
w = max(w, 0);
s = sum(w); if s > 0, w = w / s; else, w = ones(Lgrid,1)/Lgrid; end

% Reweight by isotropic-Gaussian prior if requested
if isfield(prior_cfg,'mode') && strcmpi(prior_cfg.mode, 'iso-gaussian')
    R0    = resolve_R0(prior_cfg.center);
    sigma = max(prior_cfg.sigma, 1e-12);

    prior_factor = zeros(Lgrid,1);
    for i = 1:Lgrid
        theta = so3_geodesic_angle(grid_R(:,:,i), R0);
        prior_factor(i) = exp( -0.5 * (theta/sigma)^2 );
    end
    w = w .* prior_factor;
    s = sum(w); if s > 0, w = w / s; else, w = ones(Lgrid,1)/Lgrid; end
end

end

% -------------------------- local helpers -------------------------------
function th = so3_geodesic_angle(R1, R2)
    C = R1' * R2;
    x = (trace(C) - 1) / 2;
    x = max(min(x, 1), -1);
    th = acos(x);
end

function R0 = resolve_R0(center_cfg)
    if nargin < 1 || isempty(center_cfg) || ~isfield(center_cfg,'type')
        R0 = eye(3); return;
    end
    switch lower(center_cfg.type)
        case 'identity'
            R0 = eye(3);
        case 'axis-angle'
            ax = center_cfg.axis(:).'; ax = ax ./ max(norm(ax), eps);
            R0 = axang2rotm([ax, center_cfg.angle]);
        case 'matrix'
            R0 = center_cfg.R0;
        otherwise
            error('resolve_R0: unknown center.type = %s', center_cfg.type);
    end
end
