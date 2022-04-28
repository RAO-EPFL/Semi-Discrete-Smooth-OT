function [phi_star, obj_star] = accgd(lr, param)
    % Accelerated Gradient Descent
    phi2_gd = zeros(param.N ,1);
    phi_gd = zeros(param.N, 1);
    obj_gd = zeros(param.K_agd, 1);
    t = 1;
    for i = 1 : param.K_agd
        param.eps_ = 1 / 2 / sqrt(i);
        grad_k = gradient_dual_ot(phi2_gd, param.C_X_Y, param);
        phi_gd_next = phi2_gd + lr * grad_k;
        t_next = (1 + sqrt(4 * t ^ 2 + 1)) / 2;
        phi2_gd = phi_gd_next + (t - 1) / t_next * (phi_gd_next - phi_gd);
        t = t_next;
        phi_gd = phi_gd_next;
        % evaluate objective function        
%         r_X_Y = (phi_gd - C_X_Y') / lambda;
%         r_max = max(r_X_Y, [], 1);
%         r_X_Y = r_X_Y - r_max;
%         p_star = sparsemax(r_X_Y, eta);
%         spmax = lambda + lambda * r_max + lambda * sum(r_X_Y .* p_star - p_star .^ 2 ./ repmat(eta', cont_samples, 1)', 1);
%         obj_gd(i) = nu' * phi_gd - mean(spmax);  
        obj_gd(i) = objective_comp(phi_gd, param.C_X_Y, param);
    end
    phi_star = phi_gd;
    obj_star = max(obj_gd);
end