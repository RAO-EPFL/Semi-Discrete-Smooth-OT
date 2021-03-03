%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Semi-Disrete Optimal Transport: Hardness, Regularization and Numerical
% Solution
% Authors: Bahar Taskesen, Soroosh Shafieezadeh-Abadeh, Daniel Kuhn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Suboptimality and discrepancy to $\phi^\star$ of the outputs $\bar\phi_t$
% of the Algorithm~1 for the original, entopic regularized and chi2
% divergence regularized optimal transport problems.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except method regularization_methods
close all
run params.m
fprintf(method + "-Regularized Optimal Transport Problem \n")
ind = (1:9)' * logspace(0, log10(K)-1, log10(K)); 
ind = [ind(:)', K];
sub_opt = zeros(length(ind), run_count);
phi_save = zeros(length(ind), run_count);
parfor r = 1 : run_count
    fprintf('iteration %d\n', r)
    cnt = 1;
    prm = {};
    prm.ind = ind;
    prm.K = K; % 
    prm.cont_samples = cont_samples; % Number of samples for continuous distribution
    prm.N = N; % Number of discrete samples
    prm.K_agd = K_agd; % number of iterations for accelerated gradient descent
    prm.p = p; % type of the Wasserstein distance
    prm.dist_type = dist_type; % norm type to define the cost of the Wasserstein distance
    prm.lambda = lambda; % parameter of generating function F
    prm.X = random(distX, aX, bX, d, prm.cont_samples); % samples from the continuous distribution
    prm.Y = random(distY, aY, bY, d, prm.N); % samples from the discrete distribution
    prm.C_X_Y = pdist2(prm.X', prm.Y', prm.dist_type) .^ prm.p; % cost matrix between x_i and y_j
    prm.nu = ones(prm.N, 1) / prm.N; % discrete distribution
    prm.mu = ones(prm.cont_samples, 1) / prm.cont_samples; % discretized continuous distribution
    prm.eta = ones(prm.N, 1) / prm.N;
    prm.method = method;
    if method == "chi2"
        prm.cons = cons_tikhonov;
        lr = prm.lambda * 2 * prm.N; 
        [phi_star, obj_star] = accsgd(lr, prm);
    elseif method == "entropic"
        prm.cons = cons_entropic;
        lr = prm.lambda; 
        [phi_star, obj_star] = accsgd(lr, prm);
    elseif method == "non-reg"
        prm.cons = cons_non;
        % Define variables
        phi_var = sdpvar(prm.N, 1);
        t_var = sdpvar(prm.cont_samples, 1); 
        cons = (repmat(t_var', prm.N ,1)' >= repmat(phi_var', prm.cont_samples, 1)  - prm.C_X_Y);
        % Define an objective
        %%%% The regularization term that is added here is to ensure unique
        %%%% maximizer of the dual non-regularized optimization problem
        objective = sum(prm.nu .* phi_var) - sum(t_var .* prm.mu) - 1e-5 * norm(phi_var) .^ 2;
        % Set some options for YALMIP and solver
        options = sdpsettings('verbose', 0, 'solver', 'mosek');
        % Solve the problem
        sol = optimize(cons, -objective, options);
        obj_star = value(objective); 
        phi_star = value(phi_var);
    end
    prm.obj_star = obj_star;
    prm.phi_star = phi_star;
    [sub_opt_local, phi_save_local] = sgd(prm); % Apply Algorithm~1 in the paper
    sub_opt(:, r) = sub_opt_local;
    phi_save(:, r) = vecnorm(phi_save_local - prm.phi_star, 2, 1) .^ 2;
end
save(method + "_reg_sub_opt", 'sub_opt')
save(method + "_reg_phi", 'phi_save')