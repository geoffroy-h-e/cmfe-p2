%deterministic policy function
function obj = objective_function(u, X0, B, Phi, f0_sampled, Sigma, I, Lamda)
    T = length(u);
    sum_val = 0;
    for t = 1:T
        x_t = X(t);
        sum_val = sum_val + (X0' * B * ((I - Phi)^2) * f0_sampled) - (0.5 * u' * Lamda * u);
    end
    obj = -sum_val; % Negative since we want to maximize
end