function [sys] = create_ngss()
    studentnumber = 5476747;
    a = 5;
    b = 7;
    c = 7;

    A = [0.3 + a - b, 0;1, 0.5+c];
    B = [1; 0.0];
    C = eye(2);
    D = [0; 0];
    sys = ss(A, B, C, D);
end