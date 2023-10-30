% setting key values 
T = 12; %time interval
X0 = 100000; % initial number of shares 
B = [0.3375 -0.072]; %factor loading matrix
Phi = [0.7146 0; 0 0.0353]; 
Psi = [0.0378 0; 0 0.0947]; 
Sigma = 0.0428; % var e t 
Lamda = 2.14*10^-5; % cost of transaction on a typical 1k shares
I = 0
U = 0 % do we need to initialize u? 

f0_mean = [0 0];
Omega0= [0.0412 0; 0 1.3655]; %cov of fzero 
f0_sampled = mvnrnd(f0_mean, Omega0, 1);% sample f0 to initialize code 

%deterministic policy constraints
Aeq = zeros(T, T);
beq = zeros(T, 1);

for t = 2:T
    Aeq(t, t-1:t) = [1, -1];
end

Aeq(1,1) = 1;
beq(1) = X0;
Aineq = [eye(T); -eye(T)];
bineq = zeros(2*T, 1);