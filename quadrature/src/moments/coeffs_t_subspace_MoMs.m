function  out = coeffs_t_subspace_MoMs(a_lms, b_pu, params)

rng(1)


c = .5;
r2_max = params.r2_max;
r3_max = params.r3_max;
sampling_size = params.sampling_size;
n = params.n;
L = params.L;
S = params.S;
P = params.P;
tol2 = params.tol2;
tol3 = params.tol3;
d = n^2;
n_basis = numel(a_lms);


% prepare the volume coefficients
real_t_complex = get_vol_real_t_complex(L,S);
a_lms = real_t_complex*a_lms;
A_lms_tns = vol_coeffs_vec_t_tns(a_lms, L, S);


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


% form basis matrix
grids = [kx';ky';zeros(size(ky'))]; 
grids = grids';
[r,th,phi] = Cart2Sph(grids(:,1), grids(:,2), grids(:,3));


n_vol = get_num_basis(L,S);
base_mat = zeros(n^2,n_vol);
idx = 1;

for l=0:L

    Yl = sph_harmonics(l, th, phi);
    zl = sphbes_zeros(l+1,:);
    cl = sqrt(2/c^3)./abs(sphbes(l, zl, true));

    fl = sphbes(l, r*zl/c, false)*diag(cl); fl(r>c)=0;
    for s=1:S(l+1)
        fls = fl(:,s);
        for m=-l:l
            Ylm = Yl(:,m+l+1);
            base_mat(:,idx) = (fls.*Ylm);
            idx = idx + 1;
        end
    end
end


% first moment
% SO3_rule = get_SO3_rule(L,P,1);
SO3_rule = get_exact_SO3_rule(L+2*P);
q = size(SO3_rule,1);
Psi = precompute_rot_density(SO3_rule, P);
rho = Psi*[1;b_pu];
w = SO3_rule(:,4).*rho;

M1 = 0; 
for i=1:q
    alpha = SO3_rule(i,1);
    beta = SO3_rule(i,2);
    gamma = SO3_rule(i,3);
    Rot = elr2rot(alpha, beta, gamma);
    IR = get_2D_projection(A_lms_tns, L, S, base_mat, Rot, n_basis);
    M1 = M1 + w(i)*IR;
end




% second moment
% SO3_rule = get_SO3_rule(L,P,2);
SO3_rule = get_exact_SO3_rule(2*L+2*P);
q = size(SO3_rule,1);
Psi = precompute_rot_density(SO3_rule, P);
rho = Psi*[1;b_pu];
w = SO3_rule(:,4).*rho;

Y2 = 0;
G = randn(d, sampling_size);
for i=1:q
    alpha = SO3_rule(i,1);
    beta = SO3_rule(i,2);
    gamma = SO3_rule(i,3);
    Rot = elr2rot(alpha, beta, gamma);
    IR = get_2D_projection(A_lms_tns, L, S, base_mat, Rot, n_basis);
    Y2 = Y2 + w(i)*IR*(IR'*G);
end

% obtain subspaces by SVD
[U2,S2,~] = svd(Y2,'econ'); 
s2 = diag(S2);
if sum(s2(1:r2_max).^2)/sum(s2.^2)<(1-tol2)
    r2=r2_max; 
else
    r2 = find(cumsum(s2.^2)/sum(s2.^2)>(1-tol2),1,'first');
end
U2 = U2(:,1:r2);

m2 = 0;
for i=1:q
    alpha = SO3_rule(i,1);
    beta = SO3_rule(i,2);
    gamma = SO3_rule(i,3);
    Rot = elr2rot(alpha, beta, gamma);
    IR = get_2D_projection(A_lms_tns, L, S, base_mat, Rot, n_basis);
    IR = U2'*IR;
    m2 = m2 + w(i)*(IR*IR');
end




out.m2 = m2;
out.U2 = U2;
out.s2 = s2;





% third moment 
% SO3_rule = get_SO3_rule(L,P,3);
SO3_rule = get_exact_SO3_rule(3*L+2*P);
q = size(SO3_rule,1);
Psi = precompute_rot_density(SO3_rule, P);
rho = Psi*[1;b_pu];
w = SO3_rule(:,4).*rho;

Y3 = 0;
G1 = randn(d, sampling_size);
G2 = randn(d, sampling_size);
for i=1:q
    alpha = SO3_rule(i,1);
    beta = SO3_rule(i,2);
    gamma = SO3_rule(i,3);
    Rot = elr2rot(alpha, beta, gamma);
    IR = get_2D_projection(A_lms_tns, L, S, base_mat, Rot, n_basis);
    Y3 = Y3 + w(i)*IR*((G1'*IR).*(G2'*IR)).';
end

[U3,S3,~] = svd(Y3,'econ');
s3 = diag(S3);
if sum(s3(1:r3_max).^2)/sum(s3.^2)<(1-tol3)
    r3=r3_max; 
else
    r3 = find(cumsum(s3.^2)/sum(s3.^2)>(1-tol3),1,'first');
end
U3 = U3(:,1:r3);

m3 = 0;
for i=1:q
    alpha = SO3_rule(i,1);
    beta = SO3_rule(i,2);
    gamma = SO3_rule(i,3);
    Rot = elr2rot(alpha, beta, gamma);
    IR = get_2D_projection(A_lms_tns, L, S, base_mat, Rot, n_basis);
    IR = U3'*IR;
    m3 = m3 + w(i)*tns_kron(tns_kron(IR,IR),IR);
end

out.m3 = m3;
out.U3 = U3;
out.s3 = s3;



U1 = U2;
m1 = U1'*M1; 
out.U1 = U1; 
out.m1 = m1;


end
