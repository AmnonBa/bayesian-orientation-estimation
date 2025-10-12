function [Phi_lms_nodes_MoMs, Psi_MoMs] = do_precomputations(L, S, P, bases, n, c, SO3_rule_MoMs)

% first moment
% disp('precomputing things related to M1...')
% Phi_lms_nodes = precompute_subspace_projections(bases.U1, L, S, n, c, SO3_rule_MoMs.m1);
% Phi_lms_nodes_MoMs.m1 = Phi_lms_nodes;
% Psi_MoMs.m1 = precompute_rot_density(SO3_rule_MoMs.m1, P);


% second moment 
disp('precomputing things related to M2...')
Phi_lms_nodes = precompute_subspace_projections(bases.U2, L, S, n, c, SO3_rule_MoMs.m2);
Phi_lms_nodes_MoMs.m2 = Phi_lms_nodes;
Psi_MoMs.m2 = precompute_rot_density(SO3_rule_MoMs.m2, P);


% third moment 
disp('precomputing things related to M3...')
Phi_lms_nodes = precompute_subspace_projections(bases.U3, L, S, n, c, SO3_rule_MoMs.m3);
Phi_lms_nodes_MoMs.m3 = Phi_lms_nodes;
Psi_MoMs.m3 = precompute_rot_density(SO3_rule_MoMs.m3, P);


end