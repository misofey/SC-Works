A = [0, 1, 0;
    0, 0, 1;
    -6, -4, -2];
B = [0; 0; 1];
C = [1, 0, 0];
D = 0;
sys = ss(A, B, C, D);
t = linspace(0, 1, 100).';
a = 0;
b = 0;
c = 20;
u = a.*t.*t + b.*t + c;
lsim(sys, u, t);

E = [-1/3, -5/3; -10/3, 4/3];
F = [-1, 1, -1; 2, -2, 0; 1, -1, -1];
det(F)