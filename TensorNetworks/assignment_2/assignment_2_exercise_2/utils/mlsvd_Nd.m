function [C,U]=mlsvd_Nd(T)

% MLSVD_ND computes the MLSVD of an input tensor of any order.
%
% INPUT:
%   T (N-D array): original N-Dimensional tensor 
%
% OUTPUT (in the respective order): 
%   C (N-D array): MLSVD core tensor
%   U (cell)     : MLSVD factor matrices of each mode stored in a cell,
%                  where U{1,n} gives the nth mode factor matrix
% 
% Remarks: 
%   You need to return the following variables:
U = cell(1,ndims(T)); % U{1,n} gives the nth mode factor matrix
C = zeros(size(T));


% % YOUR CODE GOES HERE


end