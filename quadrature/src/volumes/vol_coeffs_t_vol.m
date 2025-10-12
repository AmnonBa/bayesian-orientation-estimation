function V = vol_coeffs_t_vol(a_lms, n, c, L, S, is_real_coef)

if is_real_coef==true
    real_t_complex = get_vol_real_t_complex(L,S);
    a_lms = real_t_complex*a_lms;
end


if mod(n,2)==0
    x = 2*pi*((-n/2):(n/2-1)); 
else
    x = 2*pi*((-(n-1)/2):((n-1)/2)); 
end
[x,y,z] = meshgrid(x);
x = x(:); y = y(:); z = z(:);


sphbes_zeros = zeros([L+1,max(S)]);
for l=0:L
    sphbes_zeros(l+1,:) = besselzero(l+.5,max(S),1);
end


nr = ceil(1.5*n); nt = ceil(1.5*n); np = ceil(1.5*n); 
[r,th,phi,w] = spherequad(nr,nt,np,c);
[kx,ky,kz] = Sph2Cart(r,th,phi);



basis_idx = 1;
V_ft = 0;
% naively implemented, can be speeded up via vectorization
for l=0:L

    Yl = sph_harmonics(l, th, phi);
    zl = sphbes_zeros(l+1,:);
    cl = sqrt(2/c^3)./abs(sphbes(l, zl, true));

    fl = sphbes(l, r*zl/c, false)*diag(cl); 
    fl(r>c)=0;
    for s=1:S(l+1)
        fls = fl(:,s);
        for m=-l:l
            Ylm = Yl(:,m+l+1);
            V_ft = V_ft + a_lms(basis_idx)*(fls.*Ylm);
            basis_idx = basis_idx + 1;
        end
    end
end



V = finufft3d3(kx,ky,kz,w.*V_ft,1,1e-10,x,y,z);
V = real(reshape(V,[n,n,n]));

end