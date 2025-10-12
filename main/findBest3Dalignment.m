function [best_axang, vol2_rot_best, best_idx, best_corr] = findBest3Dalignment(vol1, vol2, opts)
%FINDBEST3DALIGNMENT Align vol2 to vol1 by searching random SO(3) rotations.
%
%   [best_axang, vol2_rot_best, best_idx, best_corr] = findBest3Dalignment(vol1, vol2)
%   [best_axang, vol2_rot_best, best_idx, best_corr] = findBest3Dalignment(vol1, vol2, opts)
%
% Inputs
%   vol1, vol2 : d x d x d (same size)
%   opts (optional struct):
%       .L         (int, default 3000)   number of random rotations to test
%       .interp    (char, default 'cubic')  imrotate3 interpolation
%       .mode      (char, default 'crop')   imrotate3 output mode
%       .seed      (double, default [])     RNG seed (empty = no override)
%       .normalize (logical, default true)  Pearson NCC if true; else raw dot
%
% Outputs
%   best_axang     : 1x4 [ax_x ax_y ax_z angle(rad)]
%   vol2_rot_best  : rotated vol2 at the best scoring rotation
%   best_idx       : index (1..L) of the best rotation considered
%   best_corr      : score of the best rotation (NCC or dot)

    % ---------- defaults/guards ----------
    if nargin < 3 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'L'),         opts.L = 3000;    end
    if ~isfield(opts,'interp'),    opts.interp = 'cubic'; end
    if ~isfield(opts,'mode'),      opts.mode = 'crop';    end
    if ~isfield(opts,'seed'),      opts.seed = [];        end
    if ~isfield(opts,'normalize'), opts.normalize = true; end

    if ~isequal(size(vol1), size(vol2))
        error('findBest3Dalignment:SizeMismatch','vol1 and vol2 must have the same size.');
    end
    if ~isempty(opts.seed), rng(opts.seed); end

    L = opts.L;
    axang_grid = zeros(L,4);
    scores     = -inf(L,1);

    % Precompute vol1 stats for NCC
    x = double(vol1(:));
    if opts.normalize
        x = x - mean(x);
        nx = norm(x);
        use_ncc = (nx > 0);
        if ~use_ncc
            % Fall back to raw dot if vol1 is constant after centering
            warning('findBest3Dalignment:DegenerateVol1', ...
                    'vol1 is constant after centering; falling back to raw dot product.');
        end
    else
        use_ncc = false;
    end

    % ---------- search ----------
    for l = 1:L
        % Haar draw via unit quaternion
        q = randn(1,4); q = q ./ max(norm(q), eps);
        axang = quat2axang(q);
        axang_grid(l,:) = axang;

        % Rotate vol2
        v2r = imrotate3(vol2, rad2deg(axang(4)), axang(1:3), opts.interp, opts.mode);

        % Score
        if use_ncc
            y  = double(v2r(:));
            y  = y - mean(y);
            ny = norm(y);
            if ny > 0
                s = (x' * y) / (nx * ny);
            else
                % vol2 is constant after rotation → fall back to raw dot
                s = double(vol1(:))' * double(v2r(:));
            end
        else
            s = double(vol1(:))' * double(v2r(:));
        end

        if ~isfinite(s)
            % Guard against NaN/Inf – treat as very bad score
            s = -inf;
        end

        scores(l) = s;
    end

    % ---------- pick winner & build output ----------
    [best_corr, best_idx] = max(scores);
    if ~isfinite(best_corr)
        % Extremely degenerate case: everything NaN/-Inf → pick first
        best_idx  = 1;
        best_corr = scores(1);
        warning('findBest3Dalignment:AllScoresBad', ...
                'All scores non-finite; returning the first random rotation.');
    end

    best_axang = axang_grid(best_idx,:);
    vol2_rot_best = imrotate3(vol2, rad2deg(best_axang(4)), best_axang(1:3), opts.interp, opts.mode);
end
