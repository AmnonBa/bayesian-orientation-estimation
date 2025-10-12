function rot_coef = reflect_rot_coef(rot_coef, L)


rot_coef = [1; rot_coef];

[real_2_complex,complex_2_real] = get_view_real_t_complex(L);
rot_coef=real_2_complex*rot_coef;


vec_idx = @(p, u) 2*p^2+p+u+1; 

for l=0:L
    for m=-2*l:2*l
%         sph_coef(vec_idx(l,m))=(-1)^(l-m)*sph_coef(vec_idx(l,m));
        rot_coef(vec_idx(l,m))=(-1)^(m)*rot_coef(vec_idx(l,m));
    end
end

rot_coef=real(complex_2_real*rot_coef);
rot_coef = rot_coef(2:end);



% if first_entry_remove==true 
%     rot_coef = [1; rot_coef];
% end 
% 
% 
% if is_real_coef == true 
%     [real_2_complex,complex_2_real] = get_view_real_t_complex(L);
%     rot_coef=real_2_complex*rot_coef;
% end
% 
% vec_idx = @(p, u) p^2+u+p+1; 
% 
% for l=0:L
%     for m=-l:l
% %         sph_coef(vec_idx(l,m))=(-1)^(l-m)*sph_coef(vec_idx(l,m));
%         rot_coef(vec_idx(l,m))=(-1)^(m)*rot_coef(vec_idx(l,m));
%     end
% end
% 
% if is_real_coef == true 
%     rot_coef=real(complex_2_real*rot_coef);
% end
% 
% 
% if first_entry_remove==true 
%     rot_coef = rot_coef(2:end);
% end

end