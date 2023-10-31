% setting key values 
T = 12; % time interval
X0 = 100000; % initial number of shares 

B = [0.3375; -0.072]; % factor loading matrix adjusted

Phi = [0.7146 0; 0 0.0353]; 
Psi = [0.0378 0; 0 0.0947]; 
Sigma = 0.0428; % var e t 
Lamda = 2.14*10^-5; % cost of transaction on a typical 1k shares
I = eye(2); % identity matrix

f0_mean = [0; 0];
Omega0 = [0.0412 0; 0 1.3655]; % cov of fzero 
f0_sampled = mvnrnd(f0_mean, Omega0, 1)'; % sampled f0 as column vector 

% Set the optimization options
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');

% anonymous function for optimization
objFun = @(u) objective_function(u, X0, B, Phi, f0_sampled, Lamda, I);

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

% Inequality constraints for u_t <= 0
Aineq_ut = -eye(T);
bineq_ut = zeros(T, 1);

% Explicitly set u(t) to be less than or equal to 0 for all t
Aineq_ut_explicit = -eye(T);  % Negative identity matrix
bineq_ut_explicit = zeros(T, 1);  % Vector of zeros of length T

% Inequality constraints for ensuring x(t) is less than or equal to x_0
Aineq_X_leq_X0 = tril(ones(T, T));
bineq_X_leq_X0 = zeros(T, 1);  % Since the sum of u should always be less than or equal to 0

% Inequality constraints for ensuring x(t) is always positive
Aineq_X = tril(ones(T, T));
bineq_X = repmat(X0, T, 1) - [0; cumsum(abs(u0(1:T-1)))];

% Combine all the constraints
Aineq = [Aineq_ut; Aineq_X_leq_X0; Aineq_ut_explicit];
bineq = [bineq_ut; bineq_X_leq_X0; bineq_ut_explicit];


% Call optimization solver
[u_opt, obj_val] = fmincon(objFun, u0, Aineq, bineq, Aeq, beq, [], [], [], options);

[~, X_opt] = objective_function(u_opt, X0, B, Phi, f0_sampled, Lamda, I);

disp(X_opt);
