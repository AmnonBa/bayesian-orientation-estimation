function sph_coef = rotate_sph_coef(sph_coef, Rot, L)

[real_2_complex,complex_2_real] = get_view_real_t_complex(L);
sph_coef=real_2_complex*sph_coef;


vec_idx = @(p, u) 2*p^2+p+u+1;
[alpha, beta, gamma] = rot2elr(Rot);
% WigD = cell(1,2*L+1);
% for l=0:(2*L)
%     WigD{l+1}=wignerD(l, alpha, beta, gamma); 
% end

for l=0:L
    WigD = wignerD(2*l, alpha, beta, gamma); 
    idx=vec_idx(l,-2*l):vec_idx(l,2*l);
    sph_coef(idx)=WigD.'*sph_coef(idx);
end
sph_coef=real(complex_2_real*sph_coef);


% [real_2_complex,complex_2_real] = get_view_real_t_complex(L);
% if is_real_coef == true 
%     sph_coef=real_2_complex*sph_coef;
% end
% 
% vec_idx = @(p, u) p^2+u+p+1; 
% [alpha, beta, gamma] = rot2elr(Rot);
% WigD = cell(1,L+1);
% for l=0:L
%     WigD{l+1}=wignerD(l, alpha, beta, gamma); 
% end
% 
% for l=0:L
%     idx=vec_idx(l,-l):vec_idx(l,l);
%     sph_coef(idx)=WigD{l+1}.'*sph_coef(idx);
% end
% 
% if is_real_coef == true 
%     sph_coef=real(complex_2_real*sph_coef);
% end

end