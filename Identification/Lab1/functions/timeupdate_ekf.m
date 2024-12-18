function [th_kk,P_kk] = timeupdate_ekf(th_k,P_k,Q)
    % INPUT
    % th_k   mean of filtering distribution for k|k
    % P_k    covariance matrix of filtering distribution for k|k
    % Q      process noise covariance matrix
    % OUTPUT
    % th_kk  mean of predictive distribution for k+1|k
    % P_kk   covariance of predicted distribution for k+1|k

    delta_tau = 0.5115;
     
    F_k = eye(3); 
    
    th_kk = th_k;
    th_kk(3,1) = th_k(3,1) + delta_tau;

    P_kk = F_k*P_k*F_k' + Q;

end