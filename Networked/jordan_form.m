ngss = create_ngss();


A = ngss.A;
B = ngss.B;
[V, J] = jordan(ngss.A)
Vinv = inv(V);
    
p = [-2+j, -2-j];
K = place(ngss.A, ngss.B, p);
K = [K 0];

constraint = [eye(2); eye(2)];

lb = [0;0];
ub = [1;1];


points = zeros(2, 3, 2, 2);
values = zeros(2, 3, 2);


[rows, cols, elements] = meshgrid(1:2, 1:3, 1:2);

n_matrices = 2*3*2;






matrix_array = zeros(n_matrices, 3, 3);




for row = 1:2
    for col = 1:3
        % get minimum
        fun = @(x) get_matrix(x(1), x(2), A, B, V, J, Vinv, K, row, col);
        x = fmincon(fun, [0.5, 0.5], [], [], [], [], lb, ub);
        points(row, col, :, 1) = x;
        elem = get_matrix(x(1), x(2), A, B, V, J, Vinv, K, row, col);
        values(row, col, 1) = elem;
        
        % get_maximum
        fun = @(x) -get_matrix(x(1), x(2), A, B, V, J, Vinv, K, row, col);
        x = fmincon(fun, [0.5, 0.5], [], [], [], [], lb, ub);
        points(row, col, :, 2) = x;
        element = get_matrix(x(1), x(2), A, B, V, J, Vinv, K, row, col);
        values(row, col, 2) = elem;
        
    end
end


function [elem] = get_matrix(h, tau_perc, A, B, V, J, Vinv, K, row, col)
    tau = h * tau_perc;

    expAt = @(t) expA(V, J, Vinv, t);

    Fx = expAt(h);
    Fu = (expAt(h) - expAt(h-tau)) * inv(A) * B;
    G1 = (expAt(h-tau) - eye(2)) * inv(A) * B;

    F = [Fx, Fu; zeros(1, 3)];
    G = [G1; 1];
    
    F_cl = F - G * K;
    elem = F_cl(row, col);
end

function [M] = expA(V, J, Vinv, t)
    M = V * [exp(J(1, 1) * t) 0; 0 exp(J(2, 2) * t)] * Vinv;
end
