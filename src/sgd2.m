function [sub_opt_local, phi_save_local] = sgd2(param, epochs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Semi-Disrete Optimal Transport: Hardness, Regularization and Numerical
% Solution
% Authors: Bahar Taskesen, Soroosh Shafieezadeh-Abadeh, Daniel Kuhn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementation of the Algorithm~1 in the paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ind2 = (1:9)' * logspace(0, log10(epochs*param.cont_samples)-1, log10(epochs*param.cont_samples)); 
    ind2 = [ind2(:)', epochs*param.cont_samples];
    cnt = 1;
    phi_k = zeros(param.N, 1);
    phi_ave = zeros(param.N, 1);
    phi_save_local = zeros(param.N, length(ind2));
    sub_opt_local = zeros(length(ind2), 1);
    outer_counter = 0;
    for epo = 1 : epochs 
        for k = 1 : param.cont_samples  
            x_k = param.X(:, k);
            param.eps_ = param.eps_bar / sqrt(k);
            c_x_y = pdist2(x_k', param.Y', param.dist_type) .^ param.p;
            grad_k = gradient_dual_ot(phi_k, c_x_y, param);
            phi_k = phi_k + param.cons / 4 / sqrt(k) * grad_k;
            phi_ave = (k / (k + 1)) * phi_ave + (1 / (k + 1)) * phi_k;
            outer_counter = outer_counter + 1;
            if ismember(outer_counter, ind2)
                phi_save_local(:, cnt) = phi_ave;
                integral = objective_comp(phi_ave, param.C_X_Y, param);
                sub_opt_local(cnt) = param.obj_star - integral;
                cnt = cnt + 1;
            end
        end  
    end
end
