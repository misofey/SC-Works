function [p_kf,P_kf] = KF(p_0,P_0,R,Q,y)
    % INPUT: 
    % p_0       prior mean
    % P_0       prior covariance matrix
    % R         Measurement noise covariance matrices for k = 1,...,137
    % Q         Process noise covariance matrix
    % y         Measurements for k = 1,...,137
    % OUTPUT    
    % p_kf      Mean of the filtered distribution for every k
    % P_kf      Covariance matrix of the filtered distribution for every k
    p_kf = zeros(2, 137);
    P_kf = zeros(2, 2, 137);
    p_k = p_0;
    P_k = P_0;
    
    for i = 1:137
        % prediction step
        [p_k, P_k] = timeupdate_kf(p_k, P_k, Q);
        
        % correction step
        [p_k, P_k] = measupdate_kf(p_k, P_k, R(:, :, i), y(:, i));	
        p_kf(:, i) = p_k;
        P_kf(:, :, i) = P_k;
    end
end
