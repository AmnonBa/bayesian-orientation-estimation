function vol_coef = rotate_vol_coef(vol_coef, R, L, S, is_real_coef)
if is_real_coef==true
    real_t_complex = get_vol_real_t_complex(L,S);
    vol_coef = real_t_complex*vol_coef;
end


[alpha, beta, gamma] = rot2elr(R);

idx=1;
for l=0:L
    Dl = wignerD(l,alpha,beta,gamma);
    for s=1:S(l+1)
        vol_coef(idx:(idx+2*l))=Dl'*vol_coef(idx:(idx+2*l));
        idx=idx+2*l+1;
    end
end

if is_real_coef==true
    vol_coef=real(real_t_complex\vol_coef);
end

end 



% function a_lms_rot = rotate_vol_coeffs(a_lms, R, L, S)
% 
% real_t_complex = get_vol_real_t_complex(L,S);
% [alpha, beta, gamma] = rot2elr(R);
% A_lms = vol_coeffs_vec_t_tns(real_t_complex*a_lms, L, S);
% 
% A_lms_rot = zeros(size(A_lms));
% for l=0:L
%     Dl = wignerD(l,alpha,beta,gamma);
%     A_lms_rot(l+1,1:S(l+1),1:(2*l+1)) = tns_mult(A_lms(l+1,1:S(l+1),1:(2*l+1)),3,Dl,2);
% end
% 
% a_lms_rot = real(real_t_complex\vol_coeffs_tns_t_vec(A_lms_rot, L, S));
% end
% 
% 
% 
% 
% function A_lms = vol_coeffs_vec_t_tns(a_lms, L, S)
% 
% A_lms = zeros(L+1, max(S), 2*L+1);
% 
% vec_idx = 1;
% for l=0:L
%     for s=1:S(l+1)
%         for m=-l:l
%             A_lms(l+1,s,m+l+1) = a_lms(vec_idx);
%             vec_idx = vec_idx + 1;
%         end
%     end
% end
% 
% end
% 
% 
% 
% function a_lms = vol_coeffs_tns_t_vec(A_lms, L, S)
% 
% a_lms = [];
% 
% vec_idx = 1;
% for l=0:L
%     for s=1:S(l+1)
%         for m=-l:l
%             a_lms(vec_idx) = A_lms(l+1,s,m+l+1); 
%             vec_idx = vec_idx + 1;
%         end
%     end
% end
% 
% a_lms = a_lms.';
% 
% end
% 
