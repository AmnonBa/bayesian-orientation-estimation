function rot_coef = rotate_rot_coef(rot_coef, Rot, L)


vec_idx = @(p, u) 2*p^2+p+u+1; 
[alpha,beta,gamma]=rot2elr(Rot);
rot_coef=[1;rot_coef];


[R2C,C2R] = get_view_real_t_complex(L);
rot_coef=R2C*rot_coef;

for l=0:L
    Dl = wignerD(2*l,alpha,beta,gamma);
    idx=vec_idx(l,-2*l):vec_idx(l,2*l);
%     rot_coef(idx)=conj(Dl)*rot_coef(idx);
    rot_coef(idx)=Dl.'*rot_coef(idx);
end



rot_coef=real(C2R*rot_coef);
rot_coef=rot_coef(2:end);



% vec_idx = @(p, u) p^2+u+p+1; 
% [alpha,beta,gamma]=rot2elr(Rot);
% 
% if remove_first_entry==true 
%     rot_coef=[1;rot_coef];
% end
% 
% [R2C,C2R] = get_view_real_t_complex(L);
% if is_real_coef == true
%     rot_coef=R2C*rot_coef;
% end
% 
% for l=0:L
%     Dl = wignerD(l,alpha,beta,gamma);
%     idx=vec_idx(l,-l):vec_idx(l,l);
% %     rot_coef(idx)=conj(Dl)*rot_coef(idx);
%     rot_coef(idx)=Dl.'*rot_coef(idx);
% end
% 
% 
% if is_real_coef == true
%     rot_coef=real(C2R*rot_coef);
% end
% 
% if remove_first_entry==true
%     rot_coef=rot_coef(2:end);
% end


end