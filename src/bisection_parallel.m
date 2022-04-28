function [p_bisect, delta_eps] = bisection_parallel(phi, C_X_Y, param)
    eps_ = param.eps_;
    lambda = param.lambda;
    N = param.N;
    if param.method == "hyperbolic"
        F = @(s) sinh(s./lambda - sqrt(2) + 1 + asinh(1)); % Hyperbolic entropy 
        F_i_s = @(s) min(ones(size(s)), max(zeros(size(s)), 1 - 1/N * F(-s)));
        temp = - (sqrt(2) - 1) * lambda;
        L = (exp(1/lambda) + exp(-1/lambda)) / 2 / lambda;
    elseif param.method == "entropic2" 
        F = @(s) exp(s ./ param.lambda - 1);
        F_i_s = @(s) min(ones(size(s)), max(zeros(size(s)), 1 - 1/N * F(-s)));
        temp = - lambda;
        L = 1 / N / lambda * exp(-1);
    else
        fprintf("Error : Please enter a valid regularization method ... \n");
    end
    tau_u = max(C_X_Y - phi' - temp, [], 2);
    tau_l = min(C_X_Y - phi' - temp, [], 2);
    delta_eps = eps_ / L / sqrt(N);
%     fprintf(['Number of iterations for the bisection algorithm:', ...
%         num2str(ceil(log2(max((tau_u - tau_l) / delta_eps)))),'\n'])
   
    for k = 1 : ceil(log2(max((tau_u - tau_l) / delta_eps))) 
        tau = (tau_u + tau_l) ./ 2;
        p = 1 - F_i_s(C_X_Y - phi' - tau);
        indx = find(sum(p, 2) > 1);
        tau_u(indx) = tau(indx);
        indx2 = find(sum(p, 2) <= 1);
        tau_l(indx2) = tau(indx2);
    end

    p_bisect = 1 -  F_i_s(C_X_Y - phi' - tau_l);
    p_bisect = p_bisect';
%     p_ent = exp((-C_X_Y + phi') / lambda) ./ sum(exp((-C_X_Y + phi') / lambda), 2);
%     p_ent = p_ent';
end
