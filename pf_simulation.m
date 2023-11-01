% Setting key values
T = 12; % time interval
X0 = 100000; % initial number of shares

B = [0.3375; -0.072]; % factor loading matrix adjusted
Phi = [0.7146 0; 0 0.0353];
Psi = [0.0378 0; 0 0.0947];
Sigma = 0.0428; % var e t
Lamda = 2.14e-5; % cost of transaction on a typical 1k shares
I = eye(2); % identity matrix

f0_mean = [0; 0];
Omega0 = [0.0412 0; 0 1.3655]; % cov of fzero
f0_sampled = mvnrnd(f0_mean, Omega0)'; % sampled f0 as column vector

% Number of paths for f(t)
num_paths = 100;

% Generate paths for f(t)
f_paths = zeros(2, T, num_paths);
for path_num = 1:num_paths
    f_t = f0_sampled;
    for t = 1:T
        epsilon = mvnrnd([0; 0], Psi);
        f_t = (I - Phi) * f_t + epsilon';
        f_paths(:, t, path_num) = f_t;
    end
end

% Set the optimization options
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');

% Anonymous function for optimization
objFun = @(u) pf_objective_function(u, X0, B, Lamda, f_paths);

% Set initial guess for u
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
beq(T) = -X0; % total sum to reach 0 at x(T)(on a certain tolerance level)

% Inequality constraints for u(t) <= 0
Aineq_ut = -eye(T);
bineq_ut = zeros(T, 1);
Aineq_ut2 = -eye(T);  % Negative identity matrix
bineq_ut2 = zeros(T, 1);  % Vector of zeros of length T

% Inequality constraints for ensuring x(t) is less than or equal to x_0
Aineq_X_leq_X0 = tril(ones(T, T));
bineq_X_leq_X0 = zeros(T, 1);  % Since the sum of u should always be less than or equal to 0

% Inequality constraints for ensuring x(t) is always positive
Aineq_X = tril(ones(T, T));
bineq_X = repmat(X0, T, 1) - [0; cumsum(abs(u0(1:T-1)))];

% Combine all the constraints
Aineq = [Aineq_ut; Aineq_X_leq_X0; Aineq_ut2];
bineq = [bineq_ut; bineq_X_leq_X0; bineq_ut2];

% Call optimization solver
[u_opt, obj_val] = fmincon(objFun, u0, Aineq, bineq, Aeq, beq, [], [], [], options);

% Calculate X_opt using the optimized control sequence u_opt
X_opt = pf_objective_function(u_opt, X0, B, Lamda, f_paths);

disp(u_opt);
disp(obj_val);
