%pf_simulation

% setting key values 
T = 12; 
X0 = 100000; 
B = [0.3375; -0.072];
Phi = [0.7146 0; 0 0.0353]; 
Psi = [0.0378 0; 0 0.0947];
Lamda = 2.14*10^-5; 
I = eye(2);

f0_mean = [0; 0];
Omega = [0.0412 0; 0 1.3655]; % cov of fzero 
f0_sampled = mvnrnd(f0_mean, Omega, 1)'; % sampled f0 as column vector 

% Monte Carlo
iterations_pf = 1000;
opt_u = []; % store the best control sequence
best_obj = -Inf; % initial objective value (negative infinity since we are maximizing)
obj_values = zeros(iterations_pf, 1); % Store all objective values

for i = 1:iterations_pf
    % Generate a random control sequence for this sample ensuring u(t) <= 0
    u_sample = -abs(randn(T-1, 1));
    
    % Temporary X to calculate the value of X(T-1)
    X_temp = [X0; zeros(T-1, 1)];
    for t = 1:T-1
        X_temp(t+1) = X_temp(t) + u_sample(t);
    end
    
    % Ensure X(T) = 0 by modifying the last value of u_sample
    u_sample(end) = -X_temp(end);
    
    % Evaluate the objective function for this sequence
    obj_value = pf_objective_function(u_sample, X0, B, Phi, f0_sampled, Lamda, I, Psi);
    obj_values(i) = obj_value;

    % If this sequence is better than the previous best, store it
    if obj_value > best_obj
        best_obj = obj_value;
        best_u = u_sample;
    end
end

% Compute standard error
std_dev = std(obj_values);
SE = std_dev/sqrt(iterations_pf);

% Calculate the average objective value
average_obj = sum(obj_values) / iterations_pf;

% Display results
disp('Best control sequence:');
disp(best_u);
disp(['Max objective value: ', num2str(best_obj)]);
disp(['Standard error of the max objective value: ', num2str(SE)]);
disp(['Average objective value: ', num2str(average_obj)]);
