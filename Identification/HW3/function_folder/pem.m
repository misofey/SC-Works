function [Abar,Bbar,C,D,x0, J, H] = pem(theta,A0,B0,C0,D0,x00,y,u,lambda,maxiter)
% Instructions:
% Implement your Prediction Error Method for the Output Error system here.
% Use the following function inputs and outputs.
%
% Function INPUT
% theta  Paramter vector (size: depends on your mapping choice)
% A0 Initial guess for system matrix A (matrix of size n x n)
% B0 Initial guess for system matrix B (matrix of size n x m)
% C0 Initial guess for system matrix C (matrix of size l x n)
% D0 Initial guess for system matrix D (matrix of size l x m)
% x00 Initial guess for initial state (vector of size n x one)
% u System input (matrix of size N x m)
% y System output (matrix of size N x l)
% lambda regularization parameter (scalar)
% maxiter Maximum number of iterations (scalar)
%
%
% Function OUTPUT
% Abar Estimate of system matrix A (matrix of size n x n)
% Bbar Estimate of system matrix B (matrix of size n x m)
% C Estimate of system matrix C (matrix of size l x n)
% D Estimate of system matrix D (matrix of size l x m)
% x0 Estimate of initial state (vector of size n x one)





%%%%%% YOUR CODE HERE %%%%%%
% check the given inputs for consistency
[A, B, C, D, x0] = theta2matrices(theta);
assert(isequal(A, A0))
assert(isequal(B, B0))
assert(isequal(C, C0))
assert(isequal(D, D0))
assert(isequal(x0, x00))
% traj = simsystem(A, B, C, D, x0);


p = size(theta, 1);
N = size(y, 1);
n = size(A, 1);


% Gauss-Newton

% initialize
converged = false;
maxiter_reached = false;
iter = 1;
e_hist = [];
while ~converged && ~maxiter_reached
    [A, B, C, D, x0] = theta2matrices(theta);
    [y_sim, x_sim] = simsystem(A, B, C, D, x0, u);
    E = y-y_sim;
    e_hist = cat(1, e_hist, norm(E));
    psi = setup_psi(A, C, x_sim, u);
    jaco = jacobian(psi, E);
    hess = hessian(psi);
    % size(jaco)
    % size(theta)
    % size(theta - inv(hess + lambda * eye(p)) * jaco)
    % size(inv(hess + lambda * eye(p)) * jaco)
    theta_new = theta + (hess + lambda * eye(p)) \ jaco;
    % theta_new'
    if anynan(psi)
        % psi
        disp('nans in psi')
        break
    end
    if anynan(theta_new)
        disp('nans in theta:')
        % theta_new
        % psi
        % hess
        % A
        % (hess + lambda * eye(p))
        break
    end
    % Check convergence
    if norm(theta_new - theta) < 1e-3 
        converged = true; 
    elseif iter == maxiter 
        maxiter_reached = true; 
        warning('Maximum iterations reached'); 
    end
    theta = theta_new; 
    iter = iter + 1;
end
e_hist';
theta;

%%%%%% YOUR CODE HERE %%%%%%
[Abar, Bbar, C, D, x0] = theta2matrices(theta);
J = jaco;
H = hess;
     
end

function [psi] = setup_psi(A, C, x, u)
    N = size(x, 1);
    p = 11;
    n = 4;

    psi = zeros(N, p);
    dA_dtheta = zeros(n, n, p);
    dA_dtheta(1:4, 1, 1:4) = eye(4);
    
    dB_dtheta = zeros(n, 1, p);
    dB_dtheta(2:4, 5:7) = eye(3);
    
    dx0_dtheta = zeros(n, 1, p);
    dx0_dtheta(:, 1, 8:11) = eye(4);

    for p = 1:4
        dx_dtheta = dx0_dtheta(:, 1, p);
        for j = 1:N
            psi(j, p) = C * dx_dtheta;
            p1 = A * dx_dtheta;
            p2 = dA_dtheta(:, :, p) * x(j, :)';
            dx_dtheta = p1 + p2;
            if anynan(psi)
                disp('nans in A part {i}')
            end
        end
    end
    for p = 5:7
        dx_dtheta = dx0_dtheta(:, 1, p);
        for j = 1:N
            psi(j, p) = C * dx_dtheta;
            p1 = A * dx_dtheta;
            p3 = dB_dtheta(:, :, p)*u(j);
            dx_dtheta = p1 +p3;
        end
        if anynan(psi)
            disp('nans in B part {i}')
        end
    end
    for p = 8:11
        dx_dtheta = dx0_dtheta(:, 1, p);
        for j = 1:N
            psi(j, p) = C * dx_dtheta;
            p1 = A * dx_dtheta;
            dx_dtheta = p1;
        end
        if anynan(psi)
            disp('nans in x0 part {i}')
            dx0_dtheta(:, 1, p)
            A
        end
    end
end
