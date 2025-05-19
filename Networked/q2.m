function [eigenvalues, percent_grid, h_grid] = q2(ogss, upper_limit)
    h = 0.8 * upper_limit;
    tau = 0.1 * h;

    h_min = -7;
    h_max = -0.5;
    hs = 2.^(h_min:0.125:h_max);
    nh = length(hs);

    percentages = linspace(0, 0.99, 10);
    [h_grid, percent_grid] = meshgrid(hs, percentages);

    eigenvalues = zeros(size(h_grid));
    for i = 1:length(hs)
        for j = 1:length(percentages)
            eigenvalues(j, i) = eigenvalue(percent_grid(j, i), h_grid(j, i));
        end
    end
    
    % figure;
    % contourf(percent_grid, h_grid, eigenvalues);
    % grid();
    % colorbar("Ticks", 0:0.5:6, "Direction", "reverse");

    %% sampling time selection
    h_stable = 0.55;

    p = [-2+j -2-j];
    K_original = place(ogss.A, ogss.B, p);

    ntau = 100;
    taus = linspace(0, 0.99, ntau);
    eigenvalues = zeros(size(taus));
    for i = 1:ntau
        eigenvalues(i) = eigenvalue(taus(i), h_stable);
    end

    figure;
    plot(taus*h_stable, eigenvalues);
    grid("on");
    xlabel("tau");
    ylabel("rho");
    title("Largest eigenvalue vs delay");

    max_delay_for_original_controller = max(h_stable*taus(eigenvalues < 1))

    edge_system = get_delayed_system(0.55, 0.0104);

    save("max_delay_system.mat", "edge_system");

end

function [biggest] = eigenvalue(tau_perc, h)
    tau = h * tau_perc;
    ogss = create_ogss();
    Fx = expm(ogss.A * h);
    Fu = (expm(ogss.A * h) - expm(ogss.A * (h -tau))) * inv(ogss.A) * ogss.B;
    G1 = (expm(ogss.A * (h-tau)) - eye(2)) * inv(ogss.A) * ogss.B;

    F = [Fx, Fu; zeros(1, 3)];
    G = [G1; 1];

    sys = ss(F, G, eye(3), zeros(3,1), h);

    p = [-2+j, -2-j];
    K = place(ogss.A, ogss.B, p);
    K = [K 0];
    fedback = ss(F - G*K, G, sys.C, sys.D, h);
    biggest = max(abs(pole(fedback)));
end


function [sys] = create_ogss()
    studentnumber = 5476747;
    a = 5;
    b = 7;
    c = 7;

    A = [0, 0.5 - c ; 0.2 + a - b, -1];
    B = [1.0; 0.0];
    C = eye(2);
    D = [0; 0];
    sys = ss(A, B, C, D);
end

function [sys] = get_delayed_system(tau, h)
    ogss = create_ogss();
    Fx = expm(ogss.A * h);
    Fu = (expm(ogss.A * h) - expm(ogss.A * (h -tau))) * inv(ogss.A) * ogss.B;
    G1 = (expm(ogss.A * (h-tau)) - eye(2)) * inv(ogss.A) * ogss.B;

    F = [Fx, Fu; zeros(1, 3)];
    G = [G1; 1];

    sys = ss(F, G, eye(3), zeros(3,1), h);
end