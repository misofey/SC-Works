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