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







% Gauss-Newton

% initialize
converged = false;
maxiter_reached = false;
iter = 1;

while ~converged && ~maxiter_reached 





    %%%%%% YOUR CODE HERE %%%%%%




    

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



%%%%%% YOUR CODE HERE %%%%%%



     
end