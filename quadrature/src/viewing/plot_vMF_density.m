function out = plot_vMF_density(mus,w,k)


alpha=linspace(0.1,2*pi-0.1,20)';
beta=linspace(0.1,pi-0.1,20)';
[alpha,beta]=meshgrid(alpha,beta);
alpha=alpha(:); beta=beta(:);

x=sin(beta).*cos(alpha);
y=sin(beta).*sin(alpha);
z=cos(beta);


f=get_vMF_density(mus,w,[x';y';z'],k)*4*pi;
F = scatteredInterpolant(x,y,z,f);
[X,Y,Z] = sphere(50);
C = F(X,Y,Z);

figure 
surf(X,Y,Z,C)
shading interp;
c = colorbar;
set(c,'Fontsize', 15)

out.alpha=alpha;
out.beta=beta;
out.f=f;
out.C=C;
out.X=X;
out.Y=Y;
out.Z=Z;
out.x=x;
out.y=y;
out.z=z;


F = scatteredInterpolant(x,y,z,f);
[X,Y,Z] = sphere(50);
C = F(X,Y,Z);

surf(X,Y,Z,C)
% title('Mixture of vMFs')
shading interp
c = colorbar;
set(c,'Fontsize', 15)

end