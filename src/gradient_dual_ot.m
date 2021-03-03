function grad = gradient_dual_ot(phi, C_X_Y, param)
    if param.method == "chi2"
        r_X_Y = (phi - C_X_Y') / param.lambda;
        r_max = max(r_X_Y, [], 1);
        r_X_Y = r_X_Y - r_max;
        p_star = sparsemax(r_X_Y, param.eta);   
        grad = param.nu - mean(p_star, 2);
    elseif param.method == "entropic"
        r_X_Y = (phi - C_X_Y') / param.lambda;
        r_max = max(r_X_Y, [], 1);
        r_X_Y = r_X_Y - r_max;
        p_star = exp(r_X_Y);
        p_star = p_star ./ sum(p_star, 1);
        grad = param.nu - mean(p_star, 2);
    elseif param.method == "non-reg"
        r_x_y = phi(:) - C_X_Y(:);
        grad_k = zeros(param.N, 1);
        grad_k(r_x_y == max(r_x_y)) = 1;
        grad = param.nu - grad_k - 2 * 1e-5 * phi;    
    end
end