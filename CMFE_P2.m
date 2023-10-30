% setting key values 
T = 12; % time interval
X0 = 100000; % initial number of shares 

B = [0.3375; -0.072]; % factor loading matrix adjusted to 2x1

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

% Equality and inequality constraints
Aeq = zeros(T, T);
beq = zeros(T, 1);
for t = 2:T
    Aeq(t, t-1:t) = [1, -1];
end
Aeq(1,1) = 1;
beq(1) = X0;

Aineq = [];
bineq = [];

% Call optimization solver
[u_opt, obj_val] = fmincon(objFun, u0, Aineq, bineq, Aeq, beq, [], [], [], options);
