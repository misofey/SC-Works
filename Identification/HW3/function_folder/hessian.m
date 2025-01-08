function H = hessian(psi)
% Instructions:
% This function calculates the approximated Hessian matrix in the Gauss-Newton
% algorithm
%
%  p = size(theta); 
%  N = size(y, 1); 
%  l = size(C, 1);
%
% Function INPUT
% psi derivative of estimation error wrt theta (matrix of size l*N x p)
%
% Function OUTPUT
% H Hessian matrix (matrix of size p x p)





%%%%%% YOUR CODE HERE %%%%%%
N = size(psi, 1);
H = 2/N * (psi' * psi);
end