function p_opt = sparsemax(z, eta)
    z = z';
    [K, N] = size(z);
    eta = reshape(eta, [1, N]);

    [z_sorted, sorted_indices] = sort(z, 2, 'descend'); %sorted_indices = sort(z, axis=1)[:, ::-1]
    id = sub2ind([K, N], repmat(1 : K, N, 1)', sorted_indices);

    repeated_eta = repmat(eta, K, 1);
    eta_sorted = repeated_eta(id);
    eta_cumsum = cumsum(eta_sorted, 2);

    z_eta_sorted = z_sorted .* eta_sorted;
    z_eta_cumsum = cumsum(z_eta_sorted, 2);
    k_array = 2 + eta_cumsum .* z_sorted;

    k_selected = k_array > z_eta_cumsum;
    k_max = max(k_selected .* repmat(1:N, K, 1), [], 2); 
    k_inds = sub2ind([K, N], (1:K)', k_max);
    threshold = (z_eta_cumsum(k_inds) - 2) ./ eta_cumsum(k_inds);
    p_opt = repmat(eta, K, 1) .* max(z - repmat(threshold', N, 1)', 0) / 2;
    p_opt = p_opt';
end
