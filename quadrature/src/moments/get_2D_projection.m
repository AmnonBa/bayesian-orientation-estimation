function IR = get_2D_projection(A_lms_tns, L, S, base_mat, Rot, n_basis)

[alpha, beta, gamma] = rot2elr(Rot);

A_lms_tns_rot = zeros(size(A_lms_tns));
for l=0:L
    Dl = wignerD(l,alpha,beta,gamma);
    A_lms_tns_rot(l+1,1:S(l+1),1:(2*l+1)) = tns_mult(A_lms_tns(l+1,1:S(l+1),1:(2*l+1)),3,Dl',2);
end

A_lms_rot = vol_coeffs_tns_t_vec(A_lms_tns_rot, L, S, n_basis);
IR = base_mat*A_lms_rot;


end