%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Semi-Disrete Optimal Transport: Hardness, Regularization and Numerical
% Solution
% Authors: Bahar Taskesen, Soroosh Shafieezadeh-Abadeh, Daniel Kuhn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script provides the parameters that to run wass_approx_non_reg.m,
% wass_tikhonov_regularized.m, and wass_entropic_regularized.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(112233) % seed
d = 2; % Dimension of the problem
%%% Distritbuion types for $\mu$ and $\nu$:
distX = 'Normal'; aX = 0; bX = 1;
distY = 'Uniform'; aY = -1; bY = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run_count = 20; % number of repetion
K = 1e5; % Total number of iterations for SGD
K_sgd_nnreg = 1e5;
cont_samples = 2 * K; % number of points generated from continuous distribution $\mu$
N = 5; % number of points generated from continuous distribution $\nu$
K_agd = 1e4;
epochs_sgd_opt = 5; % Total number of iterations for SGD to get the optimal solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = 1; % type of the Wasserstein distance
dist_type = 'chebychev'; % norm type to define the cost of the Wasserstein distance
cons_non = 1;
cons_entropic = 1;
cons_tikhonov = 1;
cons_hyperbolic = 1;
lambda = .1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
