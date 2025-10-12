clear 
clc 

addpath(genpath('../src/'))
addpath '../EMDB_Data/'
addpath '~/finufft/matlab'
run('~/aspire/initpath.m')

rng(1)


V = double(ReadMRC('emd_34948.map'));
V = V/norm(V(:));
n = 64;
V = vol_downsample(V,n);
pixel_size = 1.04*196/n;


L = 10;
c = 0.5;
[L,S] = get_truncate_limit(L, n, c, 1);
N_basis = get_num_basis(L,S);
P = 3;


n_vMF = 12;
k = 5;
mus = randn(3,n_vMF);
for i=1:n_vMF
    mus(:,i) = mus(:,i)/norm(mus(:,i));
end
w = ones(1,n_vMF)/n_vMF;
view_coef = vMF_t_rot_coef(mus,w,k,P);
sph_coef = rot_t_sph_coef(view_coef, P);


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



load('result_N_1e6_SNR_0.01.mat', 'out')
x3 = out.x3;
% V3 = vol_coeffs_t_vol(x3(1:N_basis), n, c, L, S, 1);
% V3 = circ_mask.*V3;
% V3 = imgaussfilt3(V3,1);
% [Rot,~,V3_aligned,reflect]=cryo_align_densities(V,V3);
load('alignment.mat','reflect', 'Rot')
% V3_1e6{1} = V3_aligned;
% res3_1e6(1) = get_fsc(V,V3_aligned,pixel_size);

view_coef_est = x3(N_basis+1:end);
if reflect==1
    view_coef_est = reflect_rot_coef(view_coef_est, P);
end
view_coef_est = rotate_rot_coef(view_coef_est, Rot, P);
sph_coef_est = rot_t_sph_coef(view_coef_est, P);



npt = 50;
alpha=linspace(0.1,2*pi-0.1,npt)';
beta=linspace(0.1,pi-0.1,npt)';
[alpha_t,beta_t]=meshgrid(alpha,beta);
alpha_t=alpha_t(:); beta_t=beta_t(:);

x=sin(beta_t).*cos(alpha_t);
y=sin(beta_t).*sin(alpha_t);
z=cos(beta_t);

f0 = get_sph_density(beta_t, alpha_t, sph_coef, P);
f = get_sph_density(beta_t, alpha_t, sph_coef_est, P);
f0 = reshape(f0,[npt,npt]);
f = reshape(f,[npt,npt]);


figure 
imagesc(beta,alpha,f)
set(gca,'fontsize',20)
xlabel('$\beta$','Interpreter','latex', 'FontSize', 30)
ylabel('$\alpha$','Interpreter','latex', 'FontSize', 30)
colorbar;
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(gcf, 'recons_viewing.pdf', '-dpdf', '-r0');
print('recons_viewing', '-depsc');



figure 
imagesc(beta,alpha,f0)
set(gca,'fontsize',20)
xlabel('$\beta$','Interpreter','latex', 'FontSize', 30)
ylabel('$\alpha$','Interpreter','latex', 'FontSize', 30)
colorbar;
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(gcf, 'sym_viewing.pdf', '-dpdf', '-r0');
print('sym_viewing', '-depsc');


figure 
imagesc(beta,alpha,abs(f-f0))
set(gca,'fontsize',20)
xlabel('$\beta$','Interpreter','latex', 'FontSize', 30)
ylabel('$\alpha$','Interpreter','latex', 'FontSize', 30)
colorbar;
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(gcf, 'err_viewing.pdf', '-dpdf', '-r0');
print('err_viewing', '-depsc');