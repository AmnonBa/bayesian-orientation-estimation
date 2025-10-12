clear 
clc 

addpath(genpath('../src/'))
addpath '../EMDB_Data/'
addpath '~/finufft/matlab'
run('~/manopt/importmanopt.m')

rng(1)


%% generate ground truth volume
V = double(ReadMRC('emd_34948.map'));
V = V/norm(V(:));
n = 64;
V = vol_downsample(V,n);

L = 12;
c = 0.5;
[L,S] = get_truncate_limit(L, n, c, 1);
N_basis = get_num_basis(L,S);
vol_coef = vol_t_vol_coeffs(V, L, S, c, 1);
V_alms = vol_coeffs_t_vol(vol_coef, n, c, L, S, 1);

%%
P = 3;
n_vMF = 12;
k = 6;
mus = randn(3,n_vMF);
for i=1:n_vMF
    mus(:,i) = mus(:,i)/norm(mus(:,i));
end
w = ones(1,n_vMF)/n_vMF;
view_coef = vMF_t_rot_coef(mus,w,k,P);


params.r2_max = 250;
params.r3_max = 130;
params.sampling_size = 300;
params.n = n;
params.L = L;
params.S = S;
params.P = P;
params.tol2 = 1e-8;
params.tol3 = 1e-6;
% out = coeffs_t_subspace_MoMs(vol_coef, view_coef, params);
% save('test3.mat', 'out')


out = run_recons(4,P);
save('test3_L=4.mat', 'out')
out = run_recons(6,P);
save('test3_L=6.mat', 'out')
out = run_recons(8,P);
save('test3_L=8.mat', 'out')
out = run_recons(10,P);
save('test3_L=10.mat', 'out')
out = run_recons(12,P);
save('test3_L=12.mat', 'out')

function out = run_recons(L,P)
    
n = 64;
c = 0.5;
[L,S] = get_truncate_limit(L, n, c, 1);
N_basis = get_num_basis(L,S);

load('test3.mat', 'out')

%% precomputation
m1 = out.m1;
m2 = out.m2;
m3 = out.m3;
SO3_rules.m2 = get_SO3_rule(13,P,1,0);
SO3_rules.m3 = get_SO3_rule(13,P,1,0);
[Phi_lms_nodes_MoMs, Psi_MoMs] = do_precomputations(L, S, P, out, n, c, SO3_rules);
l1 = 1/norm(m1(:))^2;
l2 = 1/norm(m2(:))^2;
l3 = 1/norm(m3(:))^2;

[A_view,b] = create_viewing_constrs(P);
A = sparse([zeros(size(A_view,1),N_basis) A_view]); 




%% stage 1 
options = optimoptions('fmincon','Algorithm','sqp','Display','iter', ...
'SpecifyObjectiveGradient',true,...
'StepTolerance',1e-6,'FunctionTolerance',1e-6,'StepTolerance',1e-6,...
'OptimalityTolerance',1e-6,'MaxFunctionEvaluations',1e4,'MaxIterations',1e4);


x0 = [randn(N_basis,1)*1e-6; generate_random_view_coeffs(P,10,2)];
fun = @(x) find_cost_grad(x, m1, m2, m3, l1, l2, l3, ...
    Phi_lms_nodes_MoMs, Psi_MoMs, SO3_rules);
[f,g] = fun(x0);
clear fun 


N_rep = 5;
x2_all = cell(1,N_rep);
x2cost_all = zeros(1,N_rep);
x2_time = zeros(1,N_rep);
rng(1)
for i=1:N_rep
    x0 = [randn(N_basis,1)*1e-6; generate_random_view_coeffs(P,10,2)];
    fun = @(x) find_cost_grad(x, m1, m2, m3, l1, l2, 0, ...
        Phi_lms_nodes_MoMs, Psi_MoMs, SO3_rules);
    tic 
    [x2, xcost] = fmincon(fun,x0,A,b,[],[],[],[],[],options);
    x2_time(i) = toc;
    x2cost_all(i) = xcost;
    x2_all{i} = x2;
end


out.x2_all = x2_all;
out.x2cost_all = x2cost_all;
out.x2_time = x2_time;


%% stage 2 
options = optimoptions('fmincon','Algorithm','sqp','Display','iter', ...
'SpecifyObjectiveGradient',true,...
'StepTolerance',5e-6,'FunctionTolerance',5e-6,'StepTolerance',5e-6,...
'OptimalityTolerance',5e-6,'MaxFunctionEvaluations',1e4,'MaxIterations',1e4);


[~,i] = min(x2cost_all);
x2 = x2_all{i};
fun = @(x) find_cost_grad(x, m1, m2, m3, l1, l2, l3, ...
    Phi_lms_nodes_MoMs, Psi_MoMs, SO3_rules);
tic
x3 = fmincon(fun,x2,A,b,[],[],[],[],[],options);
x3_time = toc; 



out.x2 = x2;
out.x3 = x3;
out.x3_time = x3_time;









% save('test3_L=10.mat', 'out', 'x2_all', 'x2', 'x3', 'x2_time', 'x3_time')



% run('~/aspire/initpath.m')
% pixel_size = 1.04*196/n;

% V2 = vol_coeffs_t_vol(x2(1:N_basis), n, c, L, S, 1);
% V3 = vol_coeffs_t_vol(x3(1:N_basis), n, c, L, S, 1);
% [~,~,V2_aligned,~]=cryo_align_densities(V_alms,V2);
% [~,~,V3_aligned,~]=cryo_align_densities(V_alms,V3);
% res2 = get_fsc(V_alms,V2_aligned,pixel_size);
% res3 = get_fsc(V_alms,V3_aligned,pixel_size);


% save('test3_L=10.mat', 'out', 'x2_all', 'x2', 'x3', 'x2_time', 'x3_time',...
%     'V2_aligned', 'V3_aligned', 'res2', 'res3')


% out.V2_aligned = V2_aligned;
% out.V3_aligned = V3_aligned;
% out.res2 = res2;
% out.res3 = res3;


end