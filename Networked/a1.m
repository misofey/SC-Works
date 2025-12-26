studentnumber = 5476747
a = 5;
b = 7;
c = 7;
% q1(a, b, c);
ogss = get_ogss(a, b, c);
q2(ogss, 1);
q3()
% q2(ogss, 1);

pole(edge_system)
desired_poles = [0.96, 0.95, 0.1]

K = place(edge_system.A, edge_system.B, desired_poles);
h = 0.55;

fedback = ss(edge_system.A - edge_system.B * K, edge_system.B, edge_system.C, edge_system.D, h);

pole(fedback)
step(fedback)


% h = 0.5;
% disc = c2d(ogss, h, "zoh");
% disc_fb = ss(disc.A - disc.B * K, disc.B, C, D);
% max(abs(pole(disc_fb)))
% eigenvalue(0, h a

function [sys] = get_ogss(a, b, c)
    A = [0, 0.5 - c ; 0.2 + a - b, -1];
    B = [1.0; 0.0];
    C = eye(2);
    D = [0; 0];
    sys = ss(A, B, C, D);
end
