function [p_kk,P_kk] = timeupdate_kf(p_k,P_k,Q)
    % INPUT
    % p_k   mean of filtering distribution for k|k
    % P_k   covariance matrix of filtering distribution for k|k
    % Q     process noise covariance matrix
    % OUTPUT
    % p_kk  mean of predictive distribution for k+1|k
    % P_kk  covariance of predictive distribution for k+1|k
    
    A = eye(2);

    p_kk = A*p_k;
    P_kk = A*P_k*A' + Q;

end