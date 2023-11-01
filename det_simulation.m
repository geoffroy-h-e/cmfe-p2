
% setting key values 
T = 12; % time interval
X0 = 100000; % initial number of shares 
B = [0.3375; -0.072]; % factor loading matrix adjusted
Phi = [0.7146 0; 0 0.0353]; 
Psi = [0.0378 0; 0 0.0947]; 
Sigma = 0.0428; % var e t %not useful? 
Lamda = 2.14*10^-5; % cost of transaction on a typical 1k shares
I = eye(2); % identity matrix

iterations = 100;

% Initialization
X_opt_all = zeros(T, iterations);
avg_TC_all = zeros(iterations, 1);
avg_AG_all = zeros(iterations, 1);

tic; 

for iter = 1:iterations

    f0_mean = [0; 0];
    Omega = [0.0412 0; 0 1.3655]; % cov of fzero 
    f0_sampled = mvnrnd(f0_mean, Omega, 1)'; % sampled f0 as column vector 

    % Set the optimization options, no display for better Elapsed time
    options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp'); 
    % anonymous function
    objFun = @(u) det_objective_function(u, X0, B, Phi, f0_sampled, Lamda, I);

    % Initiate U
    u0 = zeros(T, 1);

    % Equality constraints for X(t)
    Aeq = zeros(T, T);
    beq = -X0 * ones(T, 1);
    for t = 2:T
        Aeq(t, 1:t) = ones(1, t);
    end
    
    Aeq(1,1) = 1;
    beq(1) = X0;
    Aeq(T,:) = ones(1, T);
    beq(T) = -X0; % total sum to reach 0 at x(T)(on default tolerance level) 
    
    % Inequality constraints for u(t) less or equal to 0
    Aineq_ut = -eye(T);
    bineq_ut = zeros(T, 1);
    Aineq_ut2 = -eye(T);  % Negative identity matrix
    bineq_ut2 = zeros(T, 1);  % Vector of zeros
    
    % Inequality constraints for ensuring x(t) is less than or equal to x0
    Aineq_X0 = tril(ones(T, T));
    bineq_X0 = zeros(T, 1);  % Since the sum of u should always be less than or equal to 0
    
    % Inequality constraints for ensuring x(t) is always positive
    Aineq_X = tril(ones(T, T));
    bineq_X = repmat(X0, T, 1) - [0; cumsum(abs(u0(1:T-1)))];
    
    % Combine all the constraints
    Aineq = [Aineq_ut; Aineq_X0; Aineq_ut2];
    bineq = [bineq_ut; bineq_X0; bineq_ut2];

    % Call optimization solver
    [u_opt, obj_val] = fmincon(objFun, u0, Aineq, bineq, Aeq, beq, [], [], [], options);

    [~, X_opt,TC1, TC2, AG] = det_objective_function(u_opt, X0, B, Phi, f0_sampled, Lamda, I);

    % calculate statistics for Table 1
    X_opt_all(:, iter) = X_opt;
    avg_TC_all(iter) = -sum(TC2)/ T; 
    avg_AG_all(iter) = sum(AG) / T; 
end

toc

%display data 
disp('Position size:');
disp(X_opt)

disp('Trade Path:');
disp(u_opt)

disp('Mean TC over all iterations:');
disp(mean(avg_TC_all));
disp('Mean AG over all iterations:');
disp(mean(avg_AG_all));

SE_TC = std(avg_TC_all) / sqrt(length(avg_TC_all));
SE_AG = std(avg_AG_all) / sqrt(length(avg_AG_all));

disp('Standard error for Mean TC:');
disp(SE_TC);
disp('Standard error for Mean AG:');
disp(SE_AG);

