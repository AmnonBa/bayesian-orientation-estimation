function out = align_vol_coef(a_lms, a_lms_ref, L, S, is_real_coef)

real_t_complex = get_vol_real_t_complex(L,S);
a_lms=real_t_complex*a_lms;
a_lms_ref=real_t_complex*a_lms_ref;


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
problem.cost = @(Rot) norm(rotate_vol_coef(a_lms, Rot, L, S, false)-a_lms_ref)^2;

costs=zeros(1,N_rot);
for i=1:N_rot
    costs(i)=problem.cost(Rots(:,:,i));
end

[~,i]=min(costs);
Rot0 = Rots(:,:,i);
[Rot1, cost1] = trustregions(problem, Rot0);


% with reflection 
a_lms=reflect_vol_coef(a_lms, L, S, false);
problem.cost = @(Rot) norm(rotate_vol_coef(a_lms, Rot, L, S, false)-a_lms_ref)^2;

costs=zeros(1,N_rot);
for i=1:N_rot

    costs(i)=problem.cost(Rots(:,:,i));
end

[~,i]=min(costs);
Rot0 = Rots(:,:,i);
[Rot2, cost2] = trustregions(problem, Rot0);

if cost1<cost2
    cost=cost1;
    a_lms=reflect_vol_coef(a_lms, L, S, false);
    Rot=Rot1;
    reflect=false;
    a_lms_aligned=rotate_vol_coef(a_lms, Rot, L, S, false);
else
    cost=cost2;
    Rot=Rot2;
    reflect=true;
    a_lms_aligned=rotate_vol_coef(a_lms, Rot, L, S, false);
end

if is_real_coef==true
    a_lms_aligned=real(real_t_complex\a_lms_aligned);
end

out.a_lms_aligned=a_lms_aligned;
out.Rot=Rot;
out.reflect=reflect;
out.cost=cost;



end 


