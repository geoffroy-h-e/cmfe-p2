function [obj, X] = objective_function(u, X0, B, Phi, f0_sampled, Lamda, I)
    
    T = length(u);
    X = zeros(T+1, 1);
    X(1) = X0;

    sum_val = 0;
    for t = 1:T
        X(t+1) = X(t) + u(t);
        sum_val = sum_val + (X(t) * B' * ((I - Phi)^t * f0_sampled)) - (0.5 * u(t)^2 * Lamda);
    end
    obj = -sum_val; % Negative since we want to maximize
end
