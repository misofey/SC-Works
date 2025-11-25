function [B, c, output] = cpd_als_Nd(T, R, options, varargin)
% CPD_ALS_ND Computes the canonical polyadic decomposition (CPD) of a
% tensor of any order using the alternative least squares algorithm (ALS).
%
% INPUT:
%   T (I_1 x I_2 x ... x I_N): N-D data tensor 
%   R (1 x 1): number of CPD components
%   options (struct, optional) : optimization options containing:
%          - th_relerr (1 x 1): relative error threshold
%          - maxiter   (1 x 1): max number of iterations
%   variable extra inputs (1 x N cell array OR N arrays, OPTIONAL): initialization for 
%           the factor matrices
%
% OUTPUT:
%   B (1 x N cell array): cell array containing the factor matrices for all modes 
%   c (vector):  vector containing the respective weights for the CPD
%               components
%   output (struct) : optimization options containing:
%          - numiter (1 x 1): the number of iterations the algorithm ran for
%          - relerr (1 x numiter): relative error achieved, defined as Frobenius norm of 
%           the residual of the decomposition OVER Frobenius norm of the 
%           original tensor.
%
% Remarks: 
%   The order in which the factor matrices is updated is fixed: mode-1, 
%   mode-2, ..., and, lastly, mode-N. 
%   
% Authors: Sofia-Eirini Kotti (s.e.kotti@tudelft.nl)  

% Get algorithm parameters
% Number of iterations
maxiter = options.maxiter;
%
% Threshold for relative error
th_relerr = options.th_relerr;
%
init = [];
%
% Check if the initialization for the factor matrices was given
if ~isempty(varargin)
    if length(varargin) == 1    % Given as cell
        init = varargin{:}; 
    else                        % Given as matrices 
        init = varargin;
    end
end

% A useful quantity
N = ndims(T);
%
% A perhaps useful function
func = @(A) A.' * A;

% Initialize the factor matrices 
if isempty(init)    % Randomly if not initialization was given
    for idxDim = 1 : N 
        B{idxDim} = randn(size(T, idxDim), R);
    end
else                % Otherwise use the given initialization
    assert(all(cellfun(@(x) size(x,2), init) == R), ...
        "The given initialization has a different number of rank-1 components than the given R.")
    for idxDim = 1 : N 
        B{idxDim} = init{idxDim};
    end
end

% ====================== YOUR CODE HERE ======================
% You need to calculate the following variables correctly (you should comment 
% the following 3 lines out)
B = cell(1, N);
c = 1;
relerr = Inf;

% Obtain the tensor unfoldings


% ALS iterations
for idxiter = 1 : maxiter

    % Loop over dimensions (do not forget to normalize the columns!)
    for idxDim = 1 : N
        



    end
    
    % Compute the current estimate of the mode-N unfolding of the tensor
    

    % Calculate the relative error between the estimate and the true
    % mode-N unfolding of the tensor
    % Update the following line
    relerr(idxiter) = Inf;  

    % Check stopping criterion on relative error
    if relerr(idxiter) < th_relerr 
        break;
    end

end

% Add the relative error and the number of iterations to output structure
output.relerr = relerr;
output.numiter = idxiter;

% Warning if maximum number of iterations was reached
if idxiter == options.maxiter
    warning(['The ALS algorithm reached the maximum number of ' num2str(options.maxiter) ' iterations.'])
end
