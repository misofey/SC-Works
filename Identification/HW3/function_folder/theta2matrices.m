function [Abar,Bbar,C,D,x0] = theta2matrices(theta)
%%
% Function INPUT
% theta     Parameter vector (vector of size: according to the realization)
%
%
% Function OUTPUT
% Abar      System matrix A (matrix of size n x n)
% Bbar      System matrix B (matrix of size n x m)
% C         System matrix C (matrix of size l x n)
% D         System matrix D (matrix of size l x m)
% x0        Initial state (vector of size n x one)





%%%%%% YOUR CODE HERE %%%%%%
assert(isequal(size(theta), [11, 1]))
theta_a = theta(1:4);
theta_b = theta(5:7);
theta_x0 = theta(8:11);

Abar=[theta_a(1) 1 0 0;
      theta_a(2) 0 1 0;
      theta_a(3) 0 0 1;
      theta_a(4) 0 0 0];
Bbar=[ 0; theta_b(1); theta_b(2); theta_b(3)];
C=[1 0 0 0];
D=[0];
x0=[theta_x0(1); theta_x0(2); theta_x0(3); theta_x0(4)];
end