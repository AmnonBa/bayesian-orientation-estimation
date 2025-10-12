clear 
clc 



% addpath(genpath('src/'))
% addpath ~/finufft/matlab

addpath(genpath('../src/'))
addpath '../EMDB_Data/'
addpath '~/finufft/matlab'




%% generate ground truth volume
% V = double(ReadMRC('emd_34948.map'));
% V = V/norm(V(:));


filename = 'out_N_1e5_SNR_1.mat';
out = run_reconstruction(filename);
save('result_N_1e5_SNR_1.mat', 'out')

filename = 'out_N_1e5_SNR_0.1.mat';
out = run_reconstruction(filename);
save('result_N_1e5_SNR_0.1.mat', 'out')

filename = 'out_N_1e5_SNR_0.01.mat';
out = run_reconstruction(filename);
save('result_N_1e5_SNR_0.01.mat', 'out')


function out = run_reconstruction(filename)
rng(1)
n = 64;
L = 10;
c = 0.5;
[L,S] = get_truncate_limit(L, n, c, 1);
N_basis = get_num_basis(L,S);
P = 3;


% form moments 
load(filename,'out')

r1 = 246;
r2 = 246;
r3 = 102;

U1 = out.U1;
U2 = out.U2;
U3 = out.U3;

m1 = out.m1;
m2 = out.m2;
m3 = out.m3;


U1 = U1(:,1:r1);
U2 = U2(:,1:r2);
U3 = U3(:,1:r3);

m1 = m1(1:r1);
m2 = m2(1:r2,1:r2);
m3 = m3(1:r3,1:r3,1:r3);

out.U1 = U1;
out.U2 = U2;
out.U3 = U3;

out.m1 = m1;
out.m2 = m2;
out.m3 = m3;


l1 = 1/norm(m1(:))^2;
l2 = 1/norm(m2(:))^2;
l3 = 1/norm(m3(:))^2;

%% precomputaton 
SO3_rules.m2 = get_SO3_rule(15,P,1,0);
SO3_rules.m3 = get_SO3_rule(15,P,1,0);
[Phi_lms_nodes_MoMs, Psi_MoMs] = do_precomputations(L, S, P, out, n, c, SO3_rules);

%% stage one 
options = optimoptions('fmincon','Algorithm','sqp','Display','iter', ...
    'SpecifyObjectiveGradient',true,...
    'StepTolerance',1e-6,'FunctionTolerance',1e-6,'StepTolerance',1e-6,...
    'OptimalityTolerance',1e-6,'MaxFunctionEvaluations',1e4,'MaxIterations',1e4);

[A_view,b] = create_viewing_constrs(P);
A = sparse([zeros(size(A_view,1),N_basis) A_view]); 
N_rep = 5; 


x2_time = zeros(1,N_rep);
x2_all = cell(1,N_rep);
xcost_M2  = zeros(1,N_rep);
for i=1:N_rep 
    x_initial = [randn(N_basis,1)*1e-6; generate_random_view_coeffs(P,10,2)];
    fun = @(x) find_cost_grad(x, m1, m2, m3, l1, l2, 0, Phi_lms_nodes_MoMs, Psi_MoMs, SO3_rules);
    tic 
    [x, xcost] = fmincon(fun,x_initial,A,b,[],[],[],[],[],options);
    x2_time(i) = toc;
    x2_all{i} = x;
    xcost_M2(i) = xcost;
end

%% stage two 
options = optimoptions('fmincon','Algorithm','sqp','Display','iter', ...
    'SpecifyObjectiveGradient',true,...
    'StepTolerance',5e-6,'FunctionTolerance',5e-6,'StepTolerance',5e-6,...
    'OptimalityTolerance',5e-6,'MaxFunctionEvaluations',2e4,'MaxIterations',2e4);

[~,i] = min(xcost_M2);
x2 = x2_all{i};
fun = @(x) find_cost_grad(x, m1, m2, m3, l1, l2, l3, Phi_lms_nodes_MoMs, Psi_MoMs, SO3_rules);
tic 
x3 = fmincon(fun,x2,A,b,[],[],[],[],[],options);
x3_time = toc; 

out.x2 = x2;
out.x2_all = x2_all;
out.x3 = x3;

out.x2_time = x2_time;
out.x3_time = x3_time;

end

