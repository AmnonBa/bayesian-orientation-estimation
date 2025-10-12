function out = align_sph_coef(sph_coef, sph_coef_ref, P, is_real_coef)

if is_real_coef==true 
    [R2C, C2R] = get_view_real_t_complex(P);
    sph_coef = R2C*sph_coef; 
    sph_coef_ref = R2C*sph_coef_ref; 
end

alpha=linspace(0,2*pi-0.1,20);
beta=linspace(0,pi,20);
gamma=linspace(0,2*pi-0.1,20);
[alpha,beta,gamma]=meshgrid(alpha,beta,gamma);
alpha=alpha(:); beta=beta(:); gamma=gamma(:);
N_rot=numel(alpha);


Rots=zeros(3,3,N_rot);
for i=1:N_rot
    Rots(:,:,i)=Rz(alpha(i))*Ry(beta(i))*Rz(gamma(i));
end

% without reflection 
problem.M = rotationsfactory(3,1);
problem.cost = @(Rot) norm(rotate_sph_coef(sph_coef, Rot, P, false)-sph_coef_ref)^2;

costs=zeros(1,N_rot);
for i=1:N_rot
    costs(i)=problem.cost(Rots(:,:,i));
end

[~,i]=min(costs);
Rot0 = Rots(:,:,i);
[Rot1, cost1] = trustregions(problem, Rot0);

% with reflection 
sph_coef=reflect_sph_coef(sph_coef, P, false);
problem.cost = @(Rot) norm(rotate_sph_coef(sph_coef, Rot, P, false)-sph_coef_ref)^2;

costs=zeros(1,N_rot);
for i=1:N_rot
    costs(i)=problem.cost(Rots(:,:,i));
end

[~,i]=min(costs);
Rot0 = Rots(:,:,i);
[Rot2, cost2] = trustregions(problem, Rot0);

if cost1<cost2
    cost=cost1;
    sph_coef=reflect_sph_coef(sph_coef, P, false);
    Rot=Rot1;
    reflect=false;
    sph_coef_aligned=rotate_sph_coef(sph_coef, Rot, P, false);
else
    cost=cost2;
    Rot=Rot2;
    reflect=true;
    sph_coef_aligned=rotate_sph_coef(sph_coef, Rot, P, false);
end

if is_real_coef==true
    sph_coef_aligned=real(C2R*sph_coef_aligned);
end

out.sph_coef_aligned=sph_coef_aligned;
out.Rot=Rot;
out.reflect=reflect;
out.cost=cost;


end