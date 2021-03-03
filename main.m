%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Semi-Disrete Optimal Transport: Hardness, Regularization and Numerical
% Solution
% Authors: Bahar Taskesen, Soroosh Shafieezadeh-Abadeh, Daniel Kuhn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the main script to obtain the results in Figure 1. Suboptimality 
% and discrepancy to $\phi^\star$ of the outputs $\bar\phi_t$
% of the Algorithm~1 for the original, entopic regularized and chi2
% divergence regularized optimal transport problems.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
rng('default') % seed
addpath('./src')

regularization_methods = ["entropic", "chi2", "non-reg"]; % regularization methods for optimal transport
for m = 1 : length(regularization_methods)
    method = regularization_methods(m);
    run smooth_ot.m
end

%% Plot the resutls
run plot_suboptimality.m
run plot_dualvar.m

