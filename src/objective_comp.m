function obj_val = objective_comp(phi_gd, C_X_Y, param)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Semi-Disrete Optimal Transport: Hardness, Regularization and Numerical
% Solution
% Authors: Bahar Taskesen, Soroosh Shafieezadeh-Abadeh, Daniel Kuhn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of the objective of the dual smooth optimal transport problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = sqrt(2) - 1 - asinh(1);
f_hyperbolic = @(s) param.lambda * (s .* asinh(s) - sqrt(s.^2 +1) + 1 + k .* s);
f_entropic = @(s) param.lambda .* s .* log(s);
    if param.method == "chi2"
        r_X_Y = (phi_gd - C_X_Y') / param.lambda;
        r_max = max(r_X_Y, [], 1);
        r_X_Y = r_X_Y - r_max;
        p_star = sparsemax(r_X_Y, param.eta);
        spmax = param.lambda + param.lambda * r_max + param.lambda * ...
            sum(r_X_Y .* p_star - p_star .^ 2 ./ repmat(param.eta', ...
            param.cont_samples, 1)', 1);
        obj_val = param.nu' * phi_gd - mean(spmax);  
    elseif param.method == "entropic"
        r_X_Y = (phi_gd - C_X_Y') / param.lambda;
        r_max = max(r_X_Y, [], 1);
        r_X_Y = r_X_Y - r_max;
        lse = param.lambda * r_max + param.lambda * ...
            log(sum(1/param.N * exp(r_X_Y), 1));
        obj_val = param.nu' * phi_gd - mean(lse); 
    elseif param.method == "entropic2"
        p_bisect = bisection_parallel(phi_gd, C_X_Y, param);
        obj_val = param.nu' * phi_gd - mean(sum((phi_gd' - C_X_Y) .* ...
            p_bisect', 2)' - sum(1/param.N * f_entropic(param.N * p_bisect)));
    elseif param.method == "hyperbolic"
        p_bisect = bisection_parallel(phi_gd, C_X_Y, param);
        obj_val = sum(param.nu .* phi_gd) - mean(sum((phi_gd' - C_X_Y) .* ...
            p_bisect', 2)' - sum(1/param.N * f_hyperbolic(param.N * p_bisect)));
    elseif param.method == "non-reg"
        r_int = phi_gd - C_X_Y';
        r_max = max(r_int, [], 1);
        obj_val = sum(param.nu .* phi_gd) - mean(r_max);
    end
end