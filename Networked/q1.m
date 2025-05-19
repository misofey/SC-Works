function [ogss, upper_limit] = q1(a, b, c)

    h = 1;

    A = [0, 0.5 - c ; 0.2 + a - b, -1];
    B = [1.0; 0.0];
    C = eye(2)
    D = [0; 0];

    p = [-2+j, -2-j];
    K = place(A, B, p)
    ogss = ss(A, B, C, D);
    fedback = ss(A-B*K, B, C, D);

    h_min = -8;
    h_max = 0;
    hs = 2.^(linspace(h_min, h_max, 200));
    % hs = 2.^(h_min:0.06125:h_max);
    n = length(hs);

    biggest_poles = zeros(n, 1);
    for i = 1:n
        h = hs(i);
        disc = c2d(ogss, h, "zoh");
        disc_fb = ss(disc.A - disc.B * K, disc.B, C, D);
        biggest_poles(i) = max(abs(pole(disc_fb)));
    end

    figure;
    plot(hs, biggest_poles);
    grid("on");
    xlabel("h");
    ylabel("rho");
    title("Largest eigenvalue vs sampling time")
    upper_limit = max(hs(biggest_poles<1))

end

% function [sys, A, B, C, D] = ogss()
%     A = [0, 0.5 - c ; 0.2 + a - b, -1];
%     B = [1.0; 0.0];
%     C = eye(2)
%     D = [0; 0];
%     sys = ss(A, B, C, D);
% end
