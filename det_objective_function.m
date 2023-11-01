function [obj, X, TC1, TC2, AG] = det_objective_function(u, X0, B, Phi, f0_sampled, Lamda, I)
    
    T = length(u);
    X = zeros(T, 1);
    X(1) = X0;

    % To store the values 
    TC1 = zeros(T, 1);
    TC2 = zeros(T, 1); % for simulation
    AG = zeros(T, 1);

    sum_val = 0;
    for t = 1:T-1
        X(t+1) = X(t) + u(t);
        
        TC1(t) = 0.5 * u(t)^2 * Lamda; % eq to 1/2 * U(t)^T * lamda * U(t) since they are all scalars
        
        TC2(t) = u(t)^2 * Lamda;

        AG(t) = X(t) * B' * ((I - Phi)^t * f0_sampled);
        
        sum_val = sum_val + AG(t) - TC1(t);
    end
    obj = -sum_val; % Negative since we want to maximize
end