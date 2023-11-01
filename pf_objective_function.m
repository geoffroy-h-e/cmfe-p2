function obj = pf_objective_function(u, X0, B, Phi, f0_sampled, Lamda, I, Psi)
    T = length(u);
    X = zeros(T, 1);
    X(1) = X0;
    
    f = f0_sampled;

    sum_val = 0;
    for t = 1:T-1
        X(t+1) = X(t) + u(t);
        
        % Updating f for each timestep
        epsilon = mvnrnd([0; 0], Psi);
        f = (I - Phi) * f + epsilon';

        sum_val = sum_val + X(t) * B' * (I - Phi) * f - (0.5 * u(t)^2 * Lamda);
    end

    obj = -sum_val; % Negative since we want to maximize
end