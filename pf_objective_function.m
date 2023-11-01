function obj = pf_objective_function(u, X0, B, Lamda, f_paths)
    T = length(u);
    num_paths = size(f_paths, 10); % Number of paths
    X = zeros(T, 1);
    X(1) = X0;

    sum_val = 0;
    for t = 1:T-1
        X(t+1) = X(t) + u(t);

        % Calculate the expected value over multiple paths of f(t)
        expected_f = zeros(size(f_paths, 1), 1);
        for path_num = 1:num_paths
            expected_f = expected_f + f_paths(:, t, path_num);
        end
        expected_f = expected_f / num_paths;

        sum_val = sum_val + (X(t) * B' * expected_f) - (0.5 * u(t)^2 * Lamda);
    end

    obj = -sum_val; % Negative since we want to maximize
end
