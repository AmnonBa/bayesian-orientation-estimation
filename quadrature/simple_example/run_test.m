clear 
clc 

% test 

addpath(genpath('../src'))
addpath ~/finufft/matlab


N = 1e4; % sample size 
SNR = 1; % signal-to-noise ratil 
n = 64; % size of downsampled images 
use_precomputed_moments = false;

rng(1)
V = double(ReadMRC('emd_34948.map'));
V = V/norm(V(:));
D = size(V,1);
V_ds = vol_downsample(V, n);
% use a smooth approximation as ground truth 
L = 5; % truncation parameters 
c = 0.5;

[L,S] = get_truncate_limit(L, n, c, 1);
N_basis = get_num_basis(L,S);
vol_coef = vol_t_vol_coeffs(V_ds, L, S, c, 1);
V_ds = vol_coeffs_t_vol(vol_coef, n, c, L, S, 1);
V = vol_upsample(V_ds,D);


%% generate viewing angles 
P = 2; % truncation parameter for viewing direction density, corresponding to 2*P in the paper
n_vMF = 12;
k = 4;
mus = randn(3,n_vMF);
for i=1:n_vMF
    mus(:,i) = mus(:,i)/norm(mus(:,i));
end
w = ones(1,n_vMF)/n_vMF;
view_coef = vMF_t_rot_coef(mus,w,k,P);
Rots = reject_sampling_wigner(N, view_coef, P, 5);



%% form estimated moments 
params.s2 = 250;
params.s3 = 120;
params.r2_max = 250 ;
params.r3_max = 120 ;
params.tol2 = 1e-8;
params.tol3 = 1e-6;
params.sigma = 0;
params.nstd = get_estimated_std(V, SNR);
params.augsize = 5;
if use_precomputed_moments==true 
    load('out.mat','out')
else
    out = form_estimated_moments(V, n, Rots, params);
end
m1 = out.m1;
m2 = out.m2;
m3 = out.m3;

%% visualize SNR 
figure 
imagesc(sum(V,3))
colormap('gray')
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(gcf, 'clean_image.pdf', '-dpdf', '-fillpage')


figure
imagesc(sum(V,3)+randn(D,D)*params.nstd)
colormap('gray')
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(gcf, 'noisy_image.pdf', '-dpdf', '-fillpage')



%% optimization parameters 
l1 = 1/norm(m1(:))^2;
l2 = 1/norm(m2(:))^2;
l3 = 1/norm(m3(:))^2;

 
SO3_rules.m2 = get_SO3_rule(L,P,2,0);
SO3_rules.m3 = get_SO3_rule(L,P,3,0);
% precomputation
[Phi_lms_nodes_MoMs, Psi_MoMs] = do_precomputations(L, S, P, out, n, c, SO3_rules);


%% stage one optimization 
options = optimoptions('fmincon','Algorithm','sqp','Display','iter', ...
    'SpecifyObjectiveGradient',true,...
    'StepTolerance',1e-6,'FunctionTolerance',1e-6,'StepTolerance',1e-6,...
    'OptimalityTolerance',1e-6,'MaxFunctionEvaluations',1e4,'MaxIterations',1e4);

[A_view,b] = create_viewing_constrs(P);
A = sparse([zeros(size(A_view,1),N_basis) A_view]); 
N_rep = 5; 


x2_all = cell(1,N_rep);
xcost_M2  = zeros(1,N_rep);
for i=1:N_rep 
    x_initial = [randn(N_basis,1)*1e-6; generate_random_view_coeffs(P,10,2)];
    fun = @(x) find_cost_grad(x, m1, m2, m3, l1, l2, 0, Phi_lms_nodes_MoMs, Psi_MoMs, SO3_rules);
    tic 
    [x, xcost] = fmincon(fun,x_initial,A,b,[],[],[],[],[],options);
    x2_all{i} = x;
    xcost_M2(i) = xcost;
end

%% stage two optimization 
options = optimoptions('fmincon','Algorithm','sqp','Display','iter', ...
    'SpecifyObjectiveGradient',true,...
    'StepTolerance',5e-6,'FunctionTolerance',5e-6,'StepTolerance',5e-6,...
    'OptimalityTolerance',5e-6,'MaxFunctionEvaluations',2e4,'MaxIterations',2e4);

[~,i] = min(xcost_M2);
x2 = x2_all{i};
fun = @(x) find_cost_grad(x, m1, m2, m3, l1, l2, l3, Phi_lms_nodes_MoMs, Psi_MoMs, SO3_rules);
x3 = fmincon(fun,x2,A,b,[],[],[],[],[],options);

%% get resolutions 
run('~/aspire/initpath.m')
circ_mask = get_circ_mask(n);
V2 = vol_coeffs_t_vol(x2(1:N_basis), n, c, L, S, 1);
V2 = circ_mask.*V2;
V2 = imgaussfilt3(V2,1);
[~,~,V2_aligned,~]=cryo_align_densities(V_ds,V2);

V3 = vol_coeffs_t_vol(x3(1:N_basis), n, c, L, S, 1);
V3 = circ_mask.*V3;
V3 = imgaussfilt3(V3,1);
[~,~,V3_aligned,~]=cryo_align_densities(V_ds,V3);

figure 
[res2, res3] = get_fsc2(V_ds,V2_aligned,V3_aligned,1.04*D/n);
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(gcf, 'FSC.pdf', '-dpdf', '-fillpage')
WriteMRC(V_ds,1,'V_ds.mrc')
WriteMRC(V2_aligned,1,'V2.mrc')
WriteMRC(V3_aligned,1,'V3.mrc')