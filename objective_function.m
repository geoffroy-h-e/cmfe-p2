function [obj, X, TC, AG] = objective_function(u, X0, B, Phi, f0_sampled, Lamda, I)
    
    T = length(u);
    X = zeros(T, 1);
    X(1) = X0;
    TC = zeros(T, 1);  % To store the values 
    AG = zeros(T, 1);

    sum_val = 0;
    for t = 1:T-1
        X(t+1) = X(t) + u(t);
        
        TC(t) = 0.5 * u(t)^2 * Lamda;
        
        AG(t) = X(t) * B' * ((I - Phi)^t * f0_sampled);
        
        sum_val = sum_val + AG(t) - TC(t);
    end
    obj = -sum_val; % Negative since we want to maximize
end