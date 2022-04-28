function [sub_opt_local, phi_save_local] = sgd(param)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Semi-Disrete Optimal Transport: Hardness, Regularization and Numerical
% Solution
% Authors: Bahar Taskesen, Soroosh Shafieezadeh-Abadeh, Daniel Kuhn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementation of the Algorithm~1 in the paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ind = param.ind;
    phi_k = zeros(param.N, 1);
    phi_ave = zeros(param.N, 1);
    phi_save_local = zeros(param.N, length(param.ind));
    sub_opt_local = zeros(length(param.ind), 1);
    for j = 1 : length(ind)
        for k = 1 : ind(j)       
            x_k = param.X(:, k);
            param.eps_ = 1 / 2 / sqrt(k);
            c_x_y = pdist2(x_k', param.Y', param.dist_type) .^ param.p;
            grad_k = gradient_dual_ot(phi_k, c_x_y, param);
            phi_k = phi_k + param.cons / 4 / sqrt(ind(j)) * grad_k;
            phi_ave = (k / (k + 1)) * phi_ave + (1 / (k + 1)) * phi_k;
        end
        phi_save_local(:, j) = phi_ave;
        integral = objective_comp(phi_ave, param.C_X_Y, param);
        sub_opt_local(j) = param.obj_star - integral;
    end
end