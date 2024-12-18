% function [th_ekf,P_ekf] = measupdate_ekf(th_k,P_k,y_k,mic_locations,R)
%     % INPUT
%     % th_k              mean of predictive distribution for k|k-1
%     % P_k               covariance of predictive distribution for k|k-1
%     % y_k               kth measurement
%     % mic_locations     microphone locations
%     % R                 measurement noise covariance matrix
%     % OUTPUT
%     % th_ekf            mean of filtering distribution for k|k
%     % p_ekf             covariance matrix of filtering distribution for k|k
% 
%     H_k = Jacobian(th_k,mic_locations);
%     h_k = f(th_k, mic_locations);
% 
%     K_k = P_k*H_k'/(H_k*P_k*H_k' + R);
% 
%     th_ekf = th_k + K_k*(y_k - h_k);
% 
%     P_ekf = P_k - K_k*H_k*P_k;
% 
% end


function [th_ekf,P_ekf] = measupdate_ekf(th_k,P_k,y_k,mic_locations,R)
    % INPUT
    % th_k              mean of predictive distribution for k|k-1
    % P_k               covariance of predictive distribution for k|k-1
    % y_k               kth measurment
    % mic_locations     microphone locations
    % R                 measurement noise covariance matrix
    % OUTPUT
    % th_ekf            mean of filtering distribution for k|k
    % p_ekf             covariance matrix of filtering distribution for k|k

    H_k = Jacobian(th_k, mic_locations);
    h_k = f(th_k, mic_locations);
    K_k = P_k*H_k'*inv(H_k*P_k*H_k' + R);
    th_ekf1 = th_k;
    th_ekf2 = K_k*(y_k - h_k);
    th_ekf = th_ekf1 + th_ekf2;
    P_ekf = P_k - K_k*H_k*P_k;
end