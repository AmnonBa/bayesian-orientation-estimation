function [mu,alpha,beta] = plot_viewing_density(P,view_coef,is_real_coef,first_entry_remove)
n = 50;
alpha = linspace(0.01,2*pi-0.01,n)';
beta = linspace(0.1,pi,n)';
[alpha,beta] = meshgrid(alpha,beta);
alpha_t = alpha(:); beta_t=beta(:);

sph_coef = rot_t_sph_coef(view_coef,P,is_real_coef,first_entry_remove);
mu = get_sph_density(beta_t,alpha_t,sph_coef,P,is_real_coef);
mu = reshape(mu,[n,n]);

fsz=25;
fszxy=45;
figure
s=pcolor(alpha,beta,mu);
s.FaceColor = 'interp';
s.EdgeColor = 'none';
colorbar
set(gca,'fontsize',fsz)
xlabel('$\alpha$',Interpreter='latex',FontSize=fszxy)
ylabel('$\beta$',Interpreter='latex',FontSize=fszxy)
% saveas(gcf,figname)

end