h_min = 0.01
h_grid = h_min:h_min:10;
feasibles = zeros(size(h_grid));

for i = 2:length(h_grid)
    h = h_grid(i);
    tau = 0.5 * h;
    [ps, tmin, xfeas, feasible] = verify_stability_2(h);
    feasibles(i) = feasible;
end

figure;
hold on;
plot(h_grid, feasibles);
grid;
ylabel("Feasibility")
xlabel("Sampling time")
ylim([-2, 2])

    
function [F, G] = get_disc_system(h)
   ogss = create_ogss();

   F = expm(ogss.A*h);
   G = (expm(ogss.A*h) - eye(size(ogss.A))) * inv(ogss.A) * ogss.B;
end


function [ps, tmin, xfeas, feasible] = verify_stability_2(h)
    ogss = create_ogss();
    
    poles = [-2+1i -2-1i];
    K1 = place(ogss.A, ogss.B, poles);
    % K1 = 0.4 * K1;

    poles = [-1 -3];
    K2 = place(ogss.A, ogss.B, poles);

    % K2 = 0.4 * K2;
    [Fh1, Gh1] = get_disc_system(h);
    [Fh2, Gh2] = get_disc_system(2*h);
    
    h1update1 = Fh1 - Gh1 * K1;
    h1update2 = Fh1 - Gh1 * K2;

    h2update1 = Fh2;
    h2update2 = Fh2 - Gh2 * K2;

    A = zeros(4, 2, 2);
    A(1, :, :) = h1update1;
    A(2, :, :) = h2update1;
    A(3, :, :) = h1update2;
    A(4, :, :) = h2update2;

    p = 0.2;
    q = 0.3;

    transitionmatrix = [
        0 q 0  (1-q) ;
        p  0  (1-p)  0 ;
        0  q  0  (1-q) ;
        p  0  (1-p)  0];

    setlmis([]);
    
    P = [lmivar(1, [2, 1]);
        lmivar(1, [2, 1]);
        lmivar(1, [2, 1]);
        lmivar(1, [2, 1]);];
    
    for i = 1:4
        for j = 1:4
            Aj = A(j, :, :);
            Aj = reshape(Aj, 2, 2);
            lmiterm([-1, 1, 1, P(i)], eye(2), eye(2));
            lmiterm([-1, 1, 1, P(j)], -transitionmatrix(i, j) * Aj', Aj);
        end
        lmiterm([-i-4, 1, 1, P(i)], eye(2), eye(2));
    end


    lmis = getlmis();
    [tmin, xfeas] = feasp(lmis);
    feasible = sign(-tmin);
    
    ps = zeros(4, 2, 2);
    if size(feasible) ~= 1
        feasible = 0;
    end
    if feasible
        for i = 1:4
            ps(i, :, :) = dec2mat(lmis, xfeas, P(i));
        end
    end
    if feasible == 0
        feasible = -1;
    end
end