function coef_sph = rot_t_sph_coef(rot_coef, P)

rot_coef=[1;rot_coef];
[real_2_complex,complex_2_real] = get_view_real_t_complex(P);
rot_coef=real_2_complex*rot_coef;


vec_idx = @(p, u) 2*p^2+p+u+1;
coef_sph = zeros(size(rot_coef));

for p=0:P
    for u=-2*p:2*p
        coef_sph(vec_idx(p,u))=((-1)^u)*rot_coef(vec_idx(p,-u))*sqrt(4*pi/(2*2*p+1));
    end
end

coef_sph=complex_2_real*coef_sph;



% if first_entry_remove==true 
%     rot_coef=[1;rot_coef];
% end
% [real_2_complex,complex_2_real] = get_view_real_t_complex(P);
% 
% if is_real_coef==true 
%     rot_coef=real_2_complex*rot_coef;
% end
% 
% vec_idx = @(p, u) p^2+u+p+1; 
% coef_sph = zeros(size(rot_coef));
% 
% for p=0:P
%     for u=-p:p
%         coef_sph(vec_idx(p,u))=((-1)^u)*rot_coef(vec_idx(p,-u))*sqrt(4*pi/(2*p+1));
%     end
% end
% 
% if is_real_coef==true 
%     coef_sph=complex_2_real*coef_sph;
% end


end