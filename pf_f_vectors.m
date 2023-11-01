% setting key values
T = 12; % time interval
f0_mean = [0; 0];
Omega = [0.0412 0; 0 1.3655]; % cov of fzero
Phi = [0.7146 0; 0 0.0353];
Psi = [0.0378 0; 0 0.0947];
num_paths = 20;

% Initialize variables
f_paths = zeros(length(f0_mean), T, num_paths); % 3D array to store multiple paths
I = eye(size(f0_mean));

for path_num = 1:num_paths
    % Initialize f(t) for this path with the sampled f0
    f_path = mvnrnd(f0_mean, Omega)';
    f_paths(:, 1, path_num) = f_path;

    for t = 1:T-1
        % Generate epsilon at each time step
        epsilon = mvnrnd([0; 0], Psi);

        % Update f(t+1) using the equation
        f_path = (I - Phi) * f_path + epsilon';
        f_paths(:, t+1, path_num) = f_path;
    end
end
