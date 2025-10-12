clear 
clc 

addpath(genpath('../src/'))
addpath '../EMDB_Data/'
addpath '~/finufft/matlab'


V = double(ReadMRC('emd_34948.map'));
V = V/norm(V(:));
n = 64;
V = vol_downsample(V,n);
pixel_size = 1.04*196/n;

L = 12;
c = 0.5;
[L,S] = get_truncate_limit(L, n, c, 1);
N_basis = get_num_basis(L,S);
vol_coef = vol_t_vol_coeffs(V, L, S, c, 1);
V_alms = vol_coeffs_t_vol(vol_coef, n, c, L, S, 1);


run('~/aspire/initpath.m')

res2_all = zeros(1,5);
res3_all = zeros(1,5);
V2_aligned_all = cell(1,5);
V3_aligned_all = cell(1,5);

%%
load('test3_L=4.mat')

L = 4;
[L,S] = get_truncate_limit(L, n, c, 1);
N_basis = get_num_basis(L,S);

x2 = out.x2; 
V2 = vol_coeffs_t_vol(x2(1:N_basis), n, c, L, S, 1);
x3 = out.x3; 
V3 = vol_coeffs_t_vol(x3(1:N_basis), n, c, L, S, 1);

[~,~,V2_aligned,~]=cryo_align_densities(V_alms,V2);
[~,~,V3_aligned,~]=cryo_align_densities(V_alms,V3);

res2_all(1) = get_fsc(V_alms,V2_aligned,pixel_size);
res3_all(1) = get_fsc(V_alms,V3_aligned,pixel_size);

V2_aligned_all{1}  = V2_aligned;
V3_aligned_all{1}  = V3_aligned;

save('test4.mat','res2_all', 'res3_all', 'V2_aligned_all', 'V3_aligned_all')

%%
load('test3_L=6.mat')

L = 6;
[L,S] = get_truncate_limit(L, n, c, 1);
N_basis = get_num_basis(L,S);

x2 = out.x2; 
V2 = vol_coeffs_t_vol(x2(1:N_basis), n, c, L, S, 1);
x3 = out.x3; 
V3 = vol_coeffs_t_vol(x3(1:N_basis), n, c, L, S, 1);
V2 = imgaussfilt3(V2,1);
V3 = imgaussfilt3(V3,1);

[~,~,V2_aligned,~]=cryo_align_densities(V_alms,V2);
[~,~,V3_aligned,~]=cryo_align_densities(V_alms,V3);

res2_all(2) = get_fsc(V_alms,V2_aligned,pixel_size);
res3_all(2) = get_fsc(V_alms,V3_aligned,pixel_size);

V2_aligned_all{2}  = V2_aligned;
V3_aligned_all{2}  = V3_aligned;

save('test4.mat','res2_all', 'res3_all', 'V2_aligned_all', 'V3_aligned_all')

%%
load('test3_L=8.mat')

L = 8;
[L,S] = get_truncate_limit(L, n, c, 1);
N_basis = get_num_basis(L,S);

x2 = out.x2; 
V2 = vol_coeffs_t_vol(x2(1:N_basis), n, c, L, S, 1);
x3 = out.x3; 
V3 = vol_coeffs_t_vol(x3(1:N_basis), n, c, L, S, 1);
V2 = imgaussfilt3(V2,1);
V3 = imgaussfilt3(V3,1);

[~,~,V2_aligned,~]=cryo_align_densities(V_alms,V2);
[~,~,V3_aligned,~]=cryo_align_densities(V_alms,V3);

res2_all(3) = get_fsc(V_alms,V2_aligned,pixel_size);
res3_all(3) = get_fsc(V_alms,V3_aligned,pixel_size);

V2_aligned_all{3}  = V2_aligned;
V3_aligned_all{3}  = V3_aligned;

save('test4.mat','res2_all', 'res3_all', 'V2_aligned_all', 'V3_aligned_all')


%%
load('test3_L=10.mat')

L = 10;
[L,S] = get_truncate_limit(L, n, c, 1);
N_basis = get_num_basis(L,S);

x2 = out.x2; 
V2 = vol_coeffs_t_vol(x2(1:N_basis), n, c, L, S, 1);
x3 = out.x3; 
V3 = vol_coeffs_t_vol(x3(1:N_basis), n, c, L, S, 1);
V2 = imgaussfilt3(V2,1);
V3 = imgaussfilt3(V3,1);

[~,~,V2_aligned,~]=cryo_align_densities(V_alms,V2);
[~,~,V3_aligned,~]=cryo_align_densities(V_alms,V3);


res2_all(4) = get_fsc(V_alms,V2_aligned,pixel_size);
res3_all(4) = get_fsc(V_alms,V3_aligned,pixel_size);

V2_aligned_all{4}  = V2_aligned;
V3_aligned_all{4}  = V3_aligned;

save('test4.mat','res2_all', 'res3_all', 'V2_aligned_all', 'V3_aligned_all')

%%
load('test3_L=12.mat')

L = 12;
[L,S] = get_truncate_limit(L, n, c, 1);
N_basis = get_num_basis(L,S);

x2 = out.x2; 
V2 = vol_coeffs_t_vol(x2(1:N_basis), n, c, L, S, 1);
x3 = out.x3; 
V3 = vol_coeffs_t_vol(x3(1:N_basis), n, c, L, S, 1);
V2 = imgaussfilt3(V2,1);
V3 = imgaussfilt3(V3,1);

[~,~,V2_aligned,~]=cryo_align_densities(V_alms,V2);
[~,~,V3_aligned,~]=cryo_align_densities(V_alms,V3);


res2_all(5) = get_fsc(V_alms,V2_aligned,pixel_size);
res3_all(5) = get_fsc(V_alms,V3_aligned,pixel_size);

V2_aligned_all{5}  = V2_aligned;
V3_aligned_all{5}  = V3_aligned;

save('test4.mat','res2_all', 'res3_all', 'V2_aligned_all', 'V3_aligned_all')