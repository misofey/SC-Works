function [] = q3()
    h = 0.55;
    tau = 0.04;
    
    get_double_delayed_system(h, tau)

    h_min = -7;
    h_max = -0.5;
    hs = 2.^(h_min:0.125:h_max);

    percentages = linspace(0, 0.99, 100);
    [h_grid, percent_grid] = meshgrid(hs, percentages);

    eigenvalues = zeros(size(h_grid));
    for i = 1:length(hs)
        for j = 1:length(percentages)
            eigenvalues(j, i) = eigenvalue(h_grid(j, i), percent_grid(j, i));
        end
    end
    

    figure;
    hold on;
    contourf(percent_grid, h_grid, eigenvalues, 30);
    contour(percent_grid, h_grid, eigenvalues, [1 1], "r--");
    clabel([], 1);
    colorbar("Direction", "reverse");
    xlabel("Delay as fraction")
    ylabel("Sampling time")
    hold off;

    %% 3.3

    h_grid = 0.0:0.125:0.7;
    eigenvalues = zeros(size(h_grid));
    benchmark_eigenvalues = zeros(size(h_grid));
    sys = get_double_delayed_system_wo_feedback(h, tau);


    ogss = create_ogss();
     for i = 1:length(h_grid)
        h = h_grid(i);
        tau = 0.5 * h;
        sys = get_double_delayed_system(h, tau);

        benchmark_eigenvalues(i) = max(eig(sys.A));
    end
    
    h = 0.3;
    Q = eye(5);
    R = 1;
    N = zeros(5, 1);
    sys = get_double_delayed_system_wo_feedback(h, 0.5*h);
    [K, S, P] = dlqr(sys.A, sys.B, Q, R, N);

    
    h_grid = 0.0:0.0025:0.75;
    eigenvalues = zeros(size(h_grid));
    for i = 1:length(h_grid)
        h = h_grid(i);
        tau = 0.5 * h;
        sys = get_double_delayed_system_wo_feedback(h, tau);
        eigenvalues(i) = max(abs(eig(sys.A - sys.B *K)));
    end

    % minimum = min(h_grid(eigenvalues < 1));
    % maximum = max(h_grid(eigenvalues < 1));
    % interval = maximum - minimum;

    % center = -4.5;
    % offset = 0.1;
    % p1 = [center-offset center+offset];
    % p2 = [center+offset*1i center-offset * 1i];
    % 
    % K1 = place(ogss.A, ogss.B, p1);
    % K2 = place(ogss.A, ogss.B, p2);
    % 
    % K = [K1 0 K2];
    % 
    % for i = 1:length(h_grid)
    %     h = h_grid(i);
    %     tau = 0.5 * h;
    %     sys = get_double_delayed_system_wo_feedback(h, tau);
    %     eigenvalues(i) = max(abs(eig(sys.A - sys.B *K)));
    % end
    % 
    % figure;
    % hold on;
    % plot(h_grid, benchmark_eigenvalues);
    % plot(h_grid, eigenvalues);
    % legend("original controller", "redesigned");
    % grid;
    % 
    % center_grid = -10:0.25:-1;
    % intervals = zeros(size(center_grid));
    % intervals2 = zeros(size(center_grid));
    % intervals3 = zeros(size(center_grid));
    % 
    % 
    % for i = 1:length(center_grid)
    %     intervals(i) = get_stable_interval(center_grid(i), 0.1);
    %     intervals2(i) = get_stable_interval(center_grid(i), 0.6);
    %     intervals3(i) = get_stable_interval(center_grid(i), 0.3);
    % end
    figure;
    hold on;
    plot(h_grid, eigenvalues);
    grid;
    ylabel("rho")
    xlabel("Sampling time")

    % plot(center_grid, intervals2);
    % plot(center_grid, intervals3);

end

function [biggest] = eigenvalue(h, tau_perc)
    tau = tau_perc * h;
    sys = get_double_delayed_system(h, tau);
    biggest = max(abs(pole(sys)));
end

function [sys] = get_double_delayed_system(h, tau)
    ogss = create_ogss();
    
    poles = [-2+1i -2-1i];
    K1 = place(ogss.A, ogss.B, poles);

    poles = [-1 -3];
    K2 = place(ogss.A, ogss.B, poles);

    Fx = expm(ogss.A * h);
    Fu1 = (expm(ogss.A * h) - expm(ogss.A * (h -tau))) * inv(ogss.A) * ogss.B;
    G1 = (expm(ogss.A * (h-tau)) - eye(2)) * inv(ogss.A) * ogss.B;

    F = [Fx, Fu1, G1;
         0 0 0 1;
         -K2 0 0];
    G = [G1; 1; 0];

    sys = ss(F - G*[K1 0 0], G, eye(4), zeros(4,1), h);
end

function [sys] = get_double_delayed_system_wo_feedback(h, tau)
    ogss = create_ogss();

    Fx = expm(ogss.A * h);
    Fu1 = (expm(ogss.A * h) - expm(ogss.A * (h -tau))) * inv(ogss.A) * ogss.B;
    G1 = (expm(ogss.A * (h-tau)) - eye(2)) * inv(ogss.A) * ogss.B;

    F = [Fx, Fu1, zeros(2);
         zeros(1, 5);
         eye(2), zeros(2, 3)];
    G = [G1; 1; zeros(2, 1)];

    sys = ss(F, G, eye(5), zeros(5, 1), h);
end

function [interval] = get_stable_interval(center, offset)
    ogss = create_ogss();
    p1 = [center-offset center+offset];
    p2 = [center+offset*1i center-offset * 1i] * 0.8;

    K1 = place(ogss.A, ogss.B, p1)
    K2 = place(ogss.A, ogss.B, p2)

    K = [K1 0 K2];



    h_grid = 0.0:0.0025:0.5;
    eigenvalues = zeros(size(h_grid));
    for i = 1:length(h_grid)
        h = h_grid(i);
        tau = 0.5 * h;
        sys = get_double_delayed_system_wo_feedback(h, tau);
        eigenvalues(i) = max(abs(eig(sys.A - sys.B *K)));
    end

    minimum = min(h_grid(eigenvalues < 1));
    maximum = max(h_grid(eigenvalues < 1));
    interval = maximum - minimum;
end