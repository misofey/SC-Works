function [th_hat, P_hat] = NLS(yk,th_hat0,R,maxiter,mic_locations)
    % INPUT: 
    % y_k               kth measurement
    % th_hat0           initial estimate
    % R                 covariance matrix of measurement noise
    % maxiter           maximum number of iterations 
    % mic_locations     microphone locations
    
    % OUTPUT    
    % th_hat            mean of kth NLS estimate
    % P_hat             covariance matrix of kth NLSestimate
    
    b = [2.93055075098328e-06 8.69556862250426e-06 2.93055075098328e-06 ...
    -4.37180521957482e-06 1.02329067216167e-05 1.29232483949712e-05...
    -5.90914331866464e-06 -2.74318767059370e-05 ];

    counter = 0;
    th_hat = th_hat0;
    dth = 10;
    P_hat = eye(3);
    W = inv(R*R);
    while (counter < maxiter) && (norm(dth) > 1e-4)
        counter = counter+1;
        % a = [123]
        F = Jacobian(th_hat, mic_locations);
        % b = [456]
        lhs = yk.' - f(th_hat, mic_locations) - b.';
        % c = [789]
        dth = (F) \ (lhs);
        th_hat = th_hat + dth;
        % d = 123
        P_hat = P_hat + inv(F'*W*F) * P_hat;
    end
    % disp("NLS converged in "+counter+" iterations");
    % size(P_hat)
    % th_hat
    % P_hat
end
