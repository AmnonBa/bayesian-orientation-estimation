function sph_coef = reflect_sph_coef(sph_coef, L)


[real_2_complex,complex_2_real] = get_view_real_t_complex(L);
sph_coef=real_2_complex*sph_coef;


vec_idx = @(p, u) 2*p^2+p+u+1; 

for l=0:L
    for m=-2*l:2*l
%         sph_coef(vec_idx(l,m))=(-1)^(l-m)*sph_coef(vec_idx(l,m));
        sph_coef(vec_idx(l,m))=(-1)^(m)*sph_coef(vec_idx(l,m));
    end
end


sph_coef=real(complex_2_real*sph_coef);



% if is_real_coef == true 
%     [real_2_complex,complex_2_real] = get_view_real_t_complex(L);
%     sph_coef=real_2_complex*sph_coef;
% end
% 
% vec_idx = @(p, u) p^2+u+p+1; 
% 
% for l=0:L
%     for m=-l:l
% %         sph_coef(vec_idx(l,m))=(-1)^(l-m)*sph_coef(vec_idx(l,m));
%         sph_coef(vec_idx(l,m))=(-1)^(m)*sph_coef(vec_idx(l,m));
%     end
% end
% 
% if is_real_coef == true 
%     sph_coef=real(complex_2_real*sph_coef);
% end

end