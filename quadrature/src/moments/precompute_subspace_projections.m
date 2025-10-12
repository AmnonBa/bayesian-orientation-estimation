function Phi_lms_nodes = precompute_subspace_projections(U, L, S, n, c, SO3_rule)

sphbes_zeros = zeros([L+1,max(S)]);
for l=0:L
    sphbes_zeros(l+1,:) = besselzero(l+.5,max(S),1);
end

if mod(n,2)==0
    k = ((-n/2):(n/2-1))/n;  
else
    k = ((-(n-1)/2):((n-1)/2))/n; 
end
[kx,ky] = meshgrid(k);
kx = kx(:); ky = ky(:); 
[r,th,phi] = Cart2Sph(kx, ky, zeros(n^2,1));


f_l_precomputed = cell(1,L+1);
Y_l_precomputed = cell(1,L+1);

for l=0:L
    zl = sphbes_zeros(l+1,:);
    cl = sqrt(2/c^3)./abs(sphbes(l, zl, true));
    f_l = sphbes(l, r*zl/c, false);
    f_l = f_l*diag(cl); 
    f_l(r>c) = 0;
    f_l_precomputed{l+1}=f_l; 
    Y_l_precomputed{l+1} = sph_harmonics(l, th, phi);
end

real_t_complex = get_vol_real_t_complex(L,S);


q = size(SO3_rule,1);
Phi_lms_nodes = cell(1,q);


for i=1:q
    Phi_lms_nodes{i} = zeros(size(U,2),size(real_t_complex,2));
    alpha = SO3_rule(i,1);
    beta = SO3_rule(i,2);
    gamma = SO3_rule(i,3);


    idx = 1;

    for l=0:L
        D_l = wignerD(l,alpha,beta,gamma);
        Y_l=Y_l_precomputed{l+1}*D_l';


        for s=1:S(l+1)
            F_ls=f_l_precomputed{l+1}(:,s).*Y_l;
            Phi_lms_nodes{i}(:,idx:(idx+2*l))=U'*F_ls;
            idx=idx+2*l+1;
        end
    end

    Phi_lms_nodes{i} = Phi_lms_nodes{i}*real_t_complex;
    
end

end