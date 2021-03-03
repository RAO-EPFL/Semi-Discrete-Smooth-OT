function [sub_opt_local, phi_save_local] = sgd(param)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Semi-Disrete Optimal Transport: Hardness, Regularization and Numerical
% Solution
% Authors: Bahar Taskesen, Soroosh Shafieezadeh-Abadeh, Daniel Kuhn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementation of the Algorithm~1 in the paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cnt = 1;
    phi_k = zeros(param.N, 1);
    phi_ave = zeros(param.N, 1);
    phi_save_local = zeros(param.N, length(param.ind));
    sub_opt_local = zeros(length(param.ind), 1);
    for k = 1 : param.K        
        x_k = param.X(:, k);
        c_x_y = pdist2(x_k', param.Y', param.dist_type) .^ param.p;
        grad_k = gradient_dual_ot(phi_k, c_x_y, param);
        phi_k = phi_k + param.cons / 4 / sqrt(param.K) * grad_k;
        phi_ave = (k / (k + 1)) * phi_ave + (1 / (k + 1)) * phi_k;
        if ismember(k, param.ind)
            phi_save_local(:, cnt) = phi_ave;
%             r_X_Y = (phi_ave - C_X_Y') / lambda;
%             r_max = max(r_X_Y, [], 1);
%             r_X_Y = r_X_Y - r_max;
%             p_star = sparsemax(r_X_Y, eta);
%             spmax = lambda + lambda * r_max + lambda * sum(r_X_Y .* p_star - p_star .^ 2 ./ repmat(eta', cont_samples, 1)', 1);
            integral = objective_comp(phi_ave, param.C_X_Y, param);
            sub_opt_local(cnt) = param.obj_star - integral;
            cnt = cnt + 1;
        end
    end
end