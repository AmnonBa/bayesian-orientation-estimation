clear 
clc 


addpath(genpath('../src/'))
addpath '../EMDB_Data/'
addpath '~/finufft/matlab'


rng(1)
n_vMF = 12;
k = 5;
M = 5;
mus = randn(3,n_vMF);
for i=1:n_vMF
    mus(:,i) = mus(:,i)/norm(mus(:,i));
end
w = ones(1,n_vMF)/n_vMF;


npt = 50;
alpha=linspace(0.1,2*pi-0.1,npt)';
beta=linspace(0.1,pi-0.1,npt)';
[alpha_t,beta_t]=meshgrid(alpha,beta);
alpha_t=alpha_t(:); beta_t=beta_t(:);

x=sin(beta_t).*cos(alpha_t);
y=sin(beta_t).*sin(alpha_t);
z=cos(beta_t);


f = get_vMF_density(mus,w,[x';y';z'],k)*4*pi;
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
% print(gcf, 'vmf.pdf', '-dpdf', '-r0');
% print(gcf, 'vmf.pdf', '-dpdf', '-fillpage')
print('vmf', '-depsc');

V = double(ReadMRC('emd_34948.map'));
n = size(V,1);




std1 = get_estimated_std(V, 1);
std2 = get_estimated_std(V, 0.1);
std3 = get_estimated_std(V, 0.01);


I = sum(V,3);


figure 
imagesc(I+std1*randn(n,n))
colormap('gray')
% colorbar 
saveas(gca,'image_highSNR.fig')
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(gcf, 'image_highSNR.pdf', '-dpdf', '-fillpage')
print('image_highSNR', '-depsc');

figure 
imagesc(I+std2*randn(n,n))
colormap('gray')
% colorbar 
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(gcf, 'image_midSNR.pdf', '-dpdf', '-fillpage')
print('image_midSNR', '-depsc');

figure 
imagesc(I+std3*randn(n,n))
colormap('gray')
% colorbar 
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(gcf, 'image_lowSNR.pdf', '-dpdf', '-fillpage')
print('image_lowSNR', '-depsc');
