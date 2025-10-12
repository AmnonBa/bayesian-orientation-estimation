function [view_coef, view_coef_mask] = enforce_symmetry(view_coef,P)
idx = 1;
view_coef_mask = ones(size(view_coef));
for p=1:P
    for u=-p:p
        if mod(p,2)==1
            view_coef(idx) = 0;
            view_coef_mask(idx) = 0;
        end
        idx = idx+1;
    end
end
end