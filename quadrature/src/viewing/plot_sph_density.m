function mu = plot_sph_density(sph_coef, P)

n = 50;
alpha=linspace(0.1,2*pi-0.1,n)';
beta=linspace(0.1,pi-0.1,n)';
[beta_t,alpha_t]=meshgrid(beta,alpha);
alpha_t=alpha_t(:); beta_t=beta_t(:);
mu = get_sph_density(beta_t, alpha_t, sph_coef, P);

mu = reshape(mu,[n,n]);
figure
imagesc(beta,alpha,mu)
set(gca,'fontsize',20)
xlabel('$\beta$','Interpreter','latex', 'FontSize', 30)
ylabel('$\alpha$','Interpreter','latex', 'FontSize', 30)
shading interp;
colorbar;



% alpha=linspace(0.1,2*pi-0.1,30)';
% beta=linspace(0.1,pi-0.1,30)';
% [alpha,beta]=meshgrid(alpha,beta);
% alpha=alpha(:); beta=beta(:);
% 
% x=sin(beta).*cos(alpha);
% y=sin(beta).*sin(alpha);
% z=cos(beta);
% 
% 
% f=get_sph_density(beta, alpha, sph_coef, P, is_real_coef);
% F = scatteredInterpolant(x,y,z,f);
% [X,Y,Z] = sphere(50);
% C = F(X,Y,Z);
% 
% figure 
% surf(X,Y,Z,C)
% shading interp;
% c = colorbar;
% set(c,'Fontsize', 15)
% 
% out.alpha=alpha;
% out.beta=beta;
% out.f=f;
% out.C=C;
% out.X=X;
% out.Y=Y;
% out.Z=Z;
% out.x=x;
% out.y=y;
% out.z=z;
end