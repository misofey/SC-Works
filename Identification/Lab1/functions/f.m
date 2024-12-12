function ftheta = f(theta, mic_locations)
    % INPUT     
    % theta             current state estimate
    % mic_locations     microphonne locations
    % OUTPUT
    % ftheta            evaluation of f at current state estimate

    c      = 343; % in [m/s]
   
    tau = repmat(theta(3), 8, 1);
    x = repmat(theta(1), 8, 1);
    y = repmat(theta(2), 8, 1);
    xm = mic_locations(:, 1);
    ym = mic_locations(:, 2);

    ftheta = tau + 1/c * sqrt((x - xm).*(x - xm) + (y - ym).*(y - ym));
end
