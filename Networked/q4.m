function [] = q4()
%Q4 Summary of this function goes here
%   Detailed explanation goes here

    ogss = create_ogss();
    
    h1 = 0.3;
    h2 = 0.6;
    tau = 0;
    
    h1sys = get_undelayed_system(h1);
    h2sys = get_undelayed_system(h2);
    
    p1 = [-2+1i -2-1i];
    p2 = [-1 -3];
    
    K1 = place(ogss.A, ogss.B, p1)
    K2 = place(ogss.A, ogss.B, p2)
    
    % K1 = [3, -5];
    % K2 = K1;
    update_matrix_1 = (h1sys.A - h1sys.B * K1);
    update_matrix_2 = (h2sys.A - h2sys.B * K2);
    
    eig(update_matrix_1)
    eig(update_matrix_2)
    x0 = [1;1];
    nt = 20;
    t = zeros(nt+1, 1);
    x = zeros(2, nt+1);
    x(:, 1) = x0;
    for i = 2:nt+1
        if (mod(i, 2) == 0)
            t(i) = t(i-1) + h1;
            x(:, i) = update_matrix_1 * x(:, i-1);
        else
            t(i) = t(i-1) + h2;
            x(:, i) = update_matrix_2 * x(:, i-1);
        end
    end



    figure;
    hold on;
    plot(t, x(1, :), "r", "DisplayName", "x(1)_{mixed}");
    plot(t, x(2, :), "g", "DisplayName", "x(2)_{mixed}");
    xlabel("t");

    x0 = [1;1];
    nt = 20;
    t = zeros(nt+1, 1);
    x = zeros(2, nt+1);
    x(:, 1) = x0;
    for i = 2:nt+1
        if (mod(i, 2) == 0)
            t(i) = t(i-1) + h1;
            x(:, i) = update_matrix_1 * x(:, i-1);
        else
            t(i) = t(i-1) + h1;
            x(:, i) = update_matrix_1 * x(:, i-1);
        end
    end

    plot(t, x(1, :), "r--", "DisplayName", "x(1)_{h1only}");
    plot(t, x(2, :), "g--", "DisplayName", "x(2)_{h1only}");
    legend();
end

function [sys] = get_undelayed_system(h)
    ogss = create_ogss();
    Fx = expm(ogss.A*h);
    G1 = (expm(ogss.A*h) - eye(2)) * inv(ogss.A) * ogss.B;

    F = Fx;
    G = G1;

    sys = ss(F, G, eye(2), zeros(2,1), h);
end