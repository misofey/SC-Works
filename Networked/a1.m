studentnumber = 5476747
a = 5;
b = 7;
c = 7;
% q1(a, b, c);
ogss = get_ogss(a, b, c);
q2(ogss, 1);
q3()

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
