function dF = Jacobian(theta,mic_locations)
    % INPUT     
    % theta             current state estimate
    % mic_locations     microphonne locations
    % OUTPUT
    % dF                evaluation of Jacobian at current state estimate

    c      = 343; % in [m/s]
    b = [2.93055075098328e-06 8.69556862250426e-06 2.93055075098328e-06 ...
    -4.37180521957482e-06 1.02329067216167e-05 1.29232483949712e-05...
    -5.90914331866464e-06 -2.74318767059370e-05 ];

    tau = repmat(theta(3), 8, 1);
    x = repmat(theta(1), 8, 1);
    y = repmat(theta(2), 8, 1);
    xm = mic_locations(:, 1);
    ym = mic_locations(:, 2);

    dF = [(x - xm)./(sqrt((x - xm).*(x - xm) + (y - ym).*(y - ym)) * c) (y - ym)./sqrt((x - xm).*(x - xm) + (y - ym).*(y - ym))/c ones(8, 1)];  % what a gorgeous function
end