function [rotation_axis, rotation_angle, grid_R, w] = generateSO3randomRotations(N, prior_cfg)
% Build an SO(3) node set by random sampling.
%
% [axis, angle, R, w] = generateSO3randomRotations(N, sampling_cfg, prior_cfg)
%
% Inputs
%   N             : number of rotations to sample
%   prior_cfg     : struct with fields (applies to *weights*, not sampling)
%                   .mode   = 'uniform' or 'iso-gaussian'
%                   .sigma  = std (radians) for geodesic Gaussian (if 'iso-gaussian')
%                   .center = struct describing mean rotation R0:
%                             .type = 'identity' | 'axis-angle' | 'matrix'
%                               if 'axis-angle': .axis(1x3), .angle (rad)
%                               if 'matrix'    : .R0 (3x3 rotation)
%
% Outputs
%   rotation_axis  : N x 3 unit vectors
%   rotation_angle : N x 1 angles in [0, pi]
%   grid_R         : 3 x 3 x N rotation matrices
%   w              : N x 1 nonnegative weights summing to 1
%
% Notes
% - The sampling_cfg controls how *nodes* are drawn.
% - The prior_cfg controls *weights* (importance reweighting for MMSE).
%   If both are non-uniform, you intentionally bias twice.

if nargin < 2 || isempty(prior_cfg)
    prior_cfg.mode = 'uniform';
end

rotation_axis  = zeros(N, 3);
rotation_angle = zeros(N, 1);
grid_R         = zeros(3,3,N);

mode = lower(string(prior_cfg.mode));

switch mode
    case "uniform"
        % Haar via random normalized quaternion
        for i = 1:N
            q = randn(1,4);
            q = q ./ max(norm(q), eps);
            axang = quat2axang(q);          % [ax_x ax_y ax_z angle], angle in [0,pi]
            rotation_axis(i,:)  = axang(1:3);
            rotation_angle(i,1) = axang(4);
            grid_R(:,:,i)       = axang2rotm(axang);
        end

    case "iso-gaussian"
        % Axis ~ uniform on S^2; angle ~ wrapped Gaussian over [0,pi]
        eps_c = max(prior_cfg.sigma, 1e-12);   % concentration/scale
        omega = linspace(0, pi, 2000);          % angle grid
        % Wrapped Gaussian-like density on SO(3) by geodesic radius
        f = (1 - cos(omega)) .* pi .* eps_c^(-3/2) .* exp(eps_c/4 - (omega/2).^2/eps_c) ...
          .* (omega - exp(-pi^2/eps_c) .* ((omega - 2*pi).*exp(pi*omega/eps_c) + (omega + 2*pi).*exp(-pi*omega/eps_c)) ./ (2*sin(omega/2)));
        f = max(f, 0);
        F = cumsum([0, f(2:end)]);              % crude CDF
        F = F ./ max(F);

        for i = 1:N
            v = randn(1,3); v = v ./ max(norm(v), eps);
            rotation_axis(i,:) = v;

            u = rand();                          % inverse-CDF by nearest
            [~, k] = min(abs(F - u));
            theta = omega(k);
            rotation_angle(i,1) = theta;

            grid_R(:,:,i) = axang2rotm([v, theta]);
        end

    otherwise
        error('generateSO3randomRotations: unknown sampling_cfg.mode = %s', sampling_cfg.mode);
end

% Base weights: uniform over sampled nodes
w = ones(N,1) / max(N,1);

end

