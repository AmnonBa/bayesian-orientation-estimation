clear 
clc 

addpath(genpath('../src/'))
addpath '../EMDB_Data/'
addpath '~/finufft/matlab'
run('~/aspire/initpath.m')
V = double(ReadMRC('emd_34948.map'));
V = V/norm(V(:));
n = 64;
V = vol_downsample(V,n);
pixel_size = 1.04*196/n;


L = 10;
c = 0.5;
[L,S] = get_truncate_limit(L, n, c, 1);
N_basis = get_num_basis(L,S);


%%
if mod(n,2)==0
    k = ((-n/2):(n/2-1))/n;  
else
    k = ((-(n-1)/2):((n-1)/2))/n; 
end
[kx,ky,kz] = meshgrid(k);
kx = kx(:); ky = ky(:); kz = kz(:); 
Rad = (k(end)-k(1))/2;
circ_mask = sqrt(kx.^2+ky.^2+kz.^2)<=0.6*Rad;
circ_mask = reshape(circ_mask,[n,n,n]);
%%



res2_1e5 = zeros(1,3);
V2_1e5 = cell(1,3);
res3_1e5 = zeros(1,3);
V3_1e5 = cell(1,3);


%%
load('result_N_1e5_SNR_1.mat', 'out')

x2 = out.x2;
V2 = vol_coeffs_t_vol(x2(1:N_basis), n, c, L, S, 1);
V2 = circ_mask.*V2;
V2 = imgaussfilt3(V2,1);
[~,~,V2_aligned,~]=cryo_align_densities(V,V2);
V2_1e5{1} = V2_aligned;
res2_1e5(1) = get_fsc(V,V2_aligned,pixel_size);

x3 = out.x3;
V3 = vol_coeffs_t_vol(x3(1:N_basis), n, c, L, S, 1);
V3 = circ_mask.*V3;
V3 = imgaussfilt3(V3,1);
[~,~,V3_aligned,~]=cryo_align_densities(V,V3);
V3_1e5{1} = V3_aligned;
res3_1e5(1) = get_fsc(V,V3_aligned,pixel_size);

%%
load('result_N_1e5_SNR_0.1.mat', 'out')

x2 = out.x2;
V2 = vol_coeffs_t_vol(x2(1:N_basis), n, c, L, S, 1);
V2 = circ_mask.*V2;
V2 = imgaussfilt3(V2,1);
[~,~,V2_aligned,~]=cryo_align_densities(V,V2);
V2_1e5{2} = V2_aligned;
res2_1e5(2) = get_fsc(V,V2_aligned,pixel_size);


x3 = out.x3;
V3 = vol_coeffs_t_vol(x3(1:N_basis), n, c, L, S, 1);
V3 = circ_mask.*V3;
V3 = imgaussfilt3(V3,1);
[~,~,V3_aligned,~]=cryo_align_densities(V,V3);
V3_1e5{2} = V3_aligned;
res3_1e5(2) = get_fsc(V,V3_aligned,pixel_size);


%%
load('result_N_1e5_SNR_0.01.mat', 'out')

x2 = out.x2;
V2 = vol_coeffs_t_vol(x2(1:N_basis), n, c, L, S, 1);
V2 = circ_mask.*V2;
V2 = imgaussfilt3(V2,1);
[~,~,V2_aligned,~]=cryo_align_densities(V,V2);
V2_1e5{3} = V2_aligned;
res2_1e5(3) = get_fsc(V,V2_aligned,pixel_size);



x3 = out.x3;
V3 = vol_coeffs_t_vol(x3(1:N_basis), n, c, L, S, 1);
V3 = circ_mask.*V3;
V3 = imgaussfilt3(V3,1);
[~,~,V3_aligned,~]=cryo_align_densities(V,V3);
V3_1e5{3} = V3_aligned;
res3_1e5(3) = get_fsc(V,V3_aligned,pixel_size);


save('result_N_1e5.mat','V3_1e5','res3_1e5','V2_1e5','res2_1e5')





res2_1e6 = zeros(1,3);
V2_1e6 = cell(1,3);
res3_1e6 = zeros(1,3);
V3_1e6 = cell(1,3);

%%
load('result_N_1e6_SNR_1.mat', 'out')

x2 = out.x2;
V2 = vol_coeffs_t_vol(x2(1:N_basis), n, c, L, S, 1);
V2 = circ_mask.*V2;
V2 = imgaussfilt3(V2,1);
[~,~,V2_aligned,~]=cryo_align_densities(V,V2);
V2_1e6{1} = V2_aligned;
res2_1e6(1) = get_fsc(V,V2_aligned,pixel_size);

x3 = out.x3;
V3 = vol_coeffs_t_vol(x3(1:N_basis), n, c, L, S, 1);
V3 = circ_mask.*V3;
V3 = imgaussfilt3(V3,1);
[~,~,V3_aligned,~]=cryo_align_densities(V,V3);
V3_1e6{1} = V3_aligned;
res3_1e6(1) = get_fsc(V,V3_aligned,pixel_size);


%%
load('result_N_1e6_SNR_0.1.mat', 'out')

x2 = out.x2;
V2 = vol_coeffs_t_vol(x2(1:N_basis), n, c, L, S, 1);
V2 = circ_mask.*V2;
V2 = imgaussfilt3(V2,1);
[~,~,V2_aligned,~]=cryo_align_densities(V,V2);
V2_1e6{2} = V2_aligned;
res2_1e6(2) = get_fsc(V,V2_aligned,pixel_size);

x3 = out.x3;
V3 = vol_coeffs_t_vol(x3(1:N_basis), n, c, L, S, 1);
V3 = circ_mask.*V3;
V3 = imgaussfilt3(V3,1);
[~,~,V3_aligned,~]=cryo_align_densities(V,V3);
V3_1e6{2} = V3_aligned;
res3_1e6(2) = get_fsc(V,V3_aligned,pixel_size);


%%
load('result_N_1e6_SNR_0.01.mat', 'out')

x2 = out.x2;
V2 = vol_coeffs_t_vol(x2(1:N_basis), n, c, L, S, 1);
V2 = circ_mask.*V2;
V2 = imgaussfilt3(V2,1);
[~,~,V2_aligned,~]=cryo_align_densities(V,V2);
V2_1e6{3} = V2_aligned;
res2_1e6(3) = get_fsc(V,V2_aligned,pixel_size);


x3 = out.x3;
V3 = vol_coeffs_t_vol(x3(1:N_basis), n, c, L, S, 1);
V3 = circ_mask.*V3;
V3 = imgaussfilt3(V3,1);
[~,~,V3_aligned,~]=cryo_align_densities(V,V3);
V3_1e6{3} = V3_aligned;
res3_1e6(3) = get_fsc(V,V3_aligned,pixel_size);


save('result_N_1e6.mat','V3_1e6','res3_1e6','V2_1e6','res2_1e6')



res2_1e4 = zeros(1,3);
V2_1e4 = cell(1,3);
res3_1e4 = zeros(1,3);
V3_1e6 = cell(1,3);

%%
load('result_N_1e4_SNR_1.mat', 'out')

x2 = out.x2;
V2 = vol_coeffs_t_vol(x2(1:N_basis), n, c, L, S, 1);
V2 = circ_mask.*V2;
V2 = imgaussfilt3(V2,1);
[~,~,V2_aligned,~]=cryo_align_densities(V,V2);
V2_1e4{1} = V2_aligned;
res2_1e4(1) = get_fsc(V,V2_aligned,pixel_size);

x3 = out.x3;
V3 = vol_coeffs_t_vol(x3(1:N_basis), n, c, L, S, 1);
V3 = circ_mask.*V3;
V3 = imgaussfilt3(V3,1);
[~,~,V3_aligned,~]=cryo_align_densities(V,V3);
V3_1e4{1} = V3_aligned;
res3_1e4(1) = get_fsc(V,V3_aligned,pixel_size);


%%
load('result_N_1e4_SNR_0.1.mat', 'out')

x2 = out.x2;
V2 = vol_coeffs_t_vol(x2(1:N_basis), n, c, L, S, 1);
V2 = circ_mask.*V2;
V2 = imgaussfilt3(V2,1);
[~,~,V2_aligned,~]=cryo_align_densities(V,V2);
V2_1e4{2} = V2_aligned;
res2_1e4(2) = get_fsc(V,V2_aligned,pixel_size);

x3 = out.x3;
V3 = vol_coeffs_t_vol(x3(1:N_basis), n, c, L, S, 1);
V3 = circ_mask.*V3;
V3 = imgaussfilt3(V3,1);
[~,~,V3_aligned,~]=cryo_align_densities(V,V3);
V3_1e4{2} = V3_aligned;
res3_1e4(2) = get_fsc(V,V3_aligned,pixel_size);


%%
load('result_N_1e4_SNR_0.01.mat', 'out')

x2 = out.x2;
V2 = vol_coeffs_t_vol(x2(1:N_basis), n, c, L, S, 1);
V2 = circ_mask.*V2;
V2 = imgaussfilt3(V2,1);
[~,~,V2_aligned,~]=cryo_align_densities(V,V2);
V2_1e4{3} = V2_aligned;
res2_1e4(3) = get_fsc(V,V2_aligned,pixel_size);


x3 = out.x3;
V3 = vol_coeffs_t_vol(x3(1:N_basis), n, c, L, S, 1);
% V3 = circ_mask.*V3;
V3 = imgaussfilt3(V3,1);
[~,~,V3_aligned,~]=cryo_align_densities(V,V3);
V3_1e4{3} = V3_aligned;
res3_1e4(3) = get_fsc(V,V3_aligned,pixel_size);


save('result_N_1e4.mat','V3_1e4','res3_1e4','V2_1e4','res2_1e4')