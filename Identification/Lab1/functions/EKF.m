function [th_ekf,P_ekf] = EKF(th_0,P_0,y,mic_locations,R,Q)
    % INPUT 
    % th_0              prior mean
    % P_0               prior covariance matrix
    % y                 measurements for k = 1,...,137
    % mic_locations     microphone locations
    % R                 measurement noise covariance matrix
    % Q                 process noise covariance matrix
        
    % OUTPUT    
    % th_ekf            mean of the filtering distribution for every k
    % P_ekf             covariance matrix of the filtering distribution for every k    
    
    th_ekf = zeros(3, 137);
    P_ekf = zeros(3, 3, 137);
    th_k = th_0;
    P_k = P_0;

    
    for i = 1:137

        % prediction step
        [th_k, P_k] = timeupdate_ekf(th_k, P_k, Q);

        % correction step
        [th_k, P_k] = measupdate_ekf(th_k, P_k, R(:, :, i), y(:, i));	

        th_ekf(:, i) = th_k;
        P_ekf(:, :, i) = P_k;
    end
    
end



