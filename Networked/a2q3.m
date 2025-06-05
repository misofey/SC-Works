function [eigenvalues, percent_grid, h_grid] = q2()
    ngss = create_ngss()

    h_min = -7;
    h_max = -0.5;
    hs = 2.^(h_min:0.06125:h_max);

    percentages = linspace(0, 0.99, 500);
    [h_grid, percent_grid] = meshgrid(hs, percentages);

    eigenvalues = zeros(size(h_grid));
    for i = 1:length(hs)
        for j = 1:length(percentages)
            eigenvalues(j, i) = eigenvalue(percent_grid(j, i), h_grid(j, i));
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
    
    %% sampling time selection
    h_stable = 0.55;

    p = [-2+1i -2-1i];
    K_original = place(ngss.A, ngss.B, p);

    % ntau = 10;
    % taus = linspace(0, 0.15, ntau);
    % eigenvalues = zeros(size(taus));
    % max_delay_for_original_controller = max(h_stable*taus(eigenvalues < 1))

    edge_system = get_delayed_system(0.55, 0.01);

    save("max_delay_system.mat", "edge_system");
    ntau = 30;
    taus = linspace(0, 0.5, ntau);
    figure;
    hold on;
    eigenvalues = zeros(size(taus));
    for i = 1:ntau
        eigenvalues(i) = eigenvalue(taus(i)/h_stable, h_stable);
    end
    max_delay_for_original_controller = max(taus(eigenvalues < 1))


    plot(taus, eigenvalues, "DisplayName", "Original controller");
    grid("on");
    xlabel("tau");
    ylabel("rho");
    % title("Largest eigenvalue vs delay")
    % max(taus(a<1))

    pole_set_1 = [0.7+0.01i 0.7 - 0.01i 0.1]
    [a, b, c] = get_delay_response_for_poles(pole_set_1, h_stable);
    plot(taus, a, "DisplayName", "First set")
    max(taus(a<1))

    pole_set_2 = [0.65+0.15i 0.65 - 0.15i 0.3]
    [a, b, c] = get_delay_response_for_poles(pole_set_2, h_stable);
    plot(taus, a, "DisplayName", "Second set")
    max(taus(a<1))

    pole_set_3 = [0.54+0.15i 0.54-0.15i 0.3]
    [a, b, c] = get_delay_response_for_poles(pole_set_3, h_stable);
    plot(taus, a, "DisplayName", "Third set")
    max(taus(a<1))

    pole_set_4 = [0.6+0.001i 0.6-0.001i 0.599]
    [a, b, c] = get_delay_response_for_poles(pole_set_4, h_stable);
    plot(taus, a, "DisplayName", "Fourth set")
    max(taus(a<1))


    pole_set_5 = [0.5+0.18i 0.5-0.18i 0.1]
    [a, b, c] = get_delay_response_for_poles(pole_set_5, h_stable);
    plot(taus, a, "DisplayName", "Fifth set")
    max(taus(a<1))

    legend();


end

function [biggest] = eigenvalue(tau_perc, h)
    tau = h * tau_perc;
    ngss = create_ngss();
    Fx = expm(ngss.A * h);
    Fu = (expm(ngss.A * h) - expm(ngss.A * (h -tau))) * inv(ngss.A) * ngss.B;
    G1 = (expm(ngss.A * (h-tau)) - eye(2)) * inv(ngss.A) * ngss.B;

    F = [Fx, Fu; zeros(1, 3)];
    G = [G1; 1];

    sys = ss(F, G, eye(3), zeros(3,1), h);

    p = [-2+j, -2-j];
    K = place(ngss.A, ngss.B, p);
    K = [K 0];
    fedback = ss(F - G*K, G, sys.C, sys.D, h);
    biggest = max(abs(pole(fedback)));
end


function [biggest, poles] = eigenvalue_pole(tau_perc, h, K)
    tau = h * tau_perc;
    ngss = create_ngss();
    Fx = expm(ngss.A * h);
    Fu = (expm(ngss.A * h) - expm(ngss.A * (h -tau))) * inv(ngss.A) * ngss.B;
    G1 = (expm(ngss.A * (h-tau)) - eye(2)) * inv(ngss.A) * ngss.B;

    F = [Fx, Fu; zeros(1, 3)];
    G = [G1; 1];

    sys = ss(F, G, eye(3), zeros(3,1), h);
    fedback = ss(F - G*K, G, sys.C, sys.D, h);
    
    poles = pole(fedback);
    biggest = max(abs(pole(fedback)));
end

function [sys] = get_delayed_system(tau, h)
    ngss = create_ngss();
    Fx = expm(ngss.A * h);
    Fu = (expm(ngss.A * h) - expm(ngss.A * (h -tau))) * inv(ngss.A) * ngss.B;
    G1 = (expm(ngss.A * (h-tau)) - eye(2)) * inv(ngss.A) * ngss.B;

    F = [Fx, Fu; zeros(1, 3)];
    G = [G1; 1];

    sys = ss(F, G, eye(3), zeros(3,1), h);
end

function [sys] = get_fedback_system(tau, h, K)
    ngss = create_ngss();
    Fx = expm(ngss.A * h);
    Fu = (expm(ngss.A * h) - expm(ngss.A * (h -tau))) * inv(ngss.A) * ngss.B;
    G1 = (expm(ngss.A * (h-tau)) - eye(2)) * inv(ngss.A) * ngss.B;

    F = [Fx, Fu; zeros(1, 3)];
    G = [G1; 1];

    sys = ss(F, G, eye(3), zeros(3,1), h);

    % p = [-2+j, -2-j];
    % K = place(ngss.A, ngss.B, p);
    sys = ss(F - G*K, G, sys.C, sys.D, h);
end

function [eigenvalues, set_poles, K] = get_delay_response_for_poles(target_poles, h)

    ntau = 30;
    taus = linspace(0, 0.5, ntau);
    eigenvalues = zeros(size(taus));
    sz = size(taus);
    set_poles = zeros(sz(1), 3);

    sys = get_delayed_system(0, h);

    K = place(sys.A, sys.B, target_poles);
    for i = 1:ntau
        [eigenvalues(i), set_poles(i, :)] = eigenvalue_pole(taus(i), h, K);
    end


end