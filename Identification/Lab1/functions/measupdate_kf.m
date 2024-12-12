function [p_kf,P_kf] = measupdate_kf(p_k,P_k,R_k,y_k)
    % INPUT
    % p_k   mean of predictive distribution for k|k-1
    % P_k   covariance of predictive distribution for k|k-1
    % R_k   kth measurement noise covariance matrix
    % y_k   kth measurment
    % OUTPUT
    % p_kf  mean of filtering distribution for k|k
    % p_kf  covariance matrix of filtering distribution for k|k

    C = eye(2);
 
    K_k = P_k*C'/(C*P_k*C' + R_k);

    p_kf = p_k + K_k*(y_k - C*p_k);

    P_kf = P_k - P_k*C'/(C*P_k*C' + R_k)*C*P_k;

end