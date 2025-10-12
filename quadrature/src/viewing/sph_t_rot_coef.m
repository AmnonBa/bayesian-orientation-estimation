function rot_coef = sph_t_rot_coef(sph_coef, P)

[real_2_complex,complex_2_real] = get_view_real_t_complex(P);
sph_coef=real_2_complex*sph_coef;


vec_idx = @(p, u) 2*p^2+p+u+1;
rot_coef = zeros(size(sph_coef));

for p=0:P
    for u=-2*p:2*p
        rot_coef(vec_idx(p,u))=(-1)^u*sph_coef(vec_idx(p,-u))/(sqrt(4*pi/(2*2*p+1)));
    end
end

rot_coef=complex_2_real*rot_coef;

% [real_2_complex,complex_2_real] = get_view_real_t_complex(P);
% 
% if is_real_coef==true 
%     sph_coef=real_2_complex*sph_coef;
% end
% 
% vec_idx = @(p, u) p^2+u+p+1; 
% rot_coef = zeros(size(sph_coef));
% 
% for p=0:P
%     for u=-p:p
%         rot_coef(vec_idx(p,u))=(-1)^u*sph_coef(vec_idx(p,-u))/(sqrt(4*pi/(2*p+1)));
%     end
% end
% 
% if is_real_coef==true 
%     rot_coef=complex_2_real*rot_coef;
% end


end