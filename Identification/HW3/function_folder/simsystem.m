function [y, x] = simsystem(A, B, C, D, x0, u)
% Instructions:
% Simulating a linear dynamic system given input u, matrices A,B,C,D ,and
% initial condition x(0)
%
 n = size(A, 1); 
 m = size(B, 2); 
 l = size(C, 1); 
 N = size(u, 1); 
%
% Function INPUT
% A system matrix (matrix of size n x n)
% B system matrix (matrix of size n x m)
% C system matrix (matrix of size l x n)
% D system matrix (matrix of size l x m)
% x0 initial state (vector of size n x one)
% u system input (matrix of size N x m)
%
% Function OUTPUT
% x state of system (matrix of size N x n)
% y system output (matrix of size N x l)




%%%%%% YOUR CODE HERE %%%%%%
y = zeros(l, N);
x = zeros(n, N+1);

x(:, 1) = x0;

for i = 1:size(u, 1)
    x(:, i+1) = (A * x(:, i) + B * u(i, :));
    y(:, i) = (C * x(:, i) + D * u(i, :));
end

x = x(:, 1:N)';
y = y';
end