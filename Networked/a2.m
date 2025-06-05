h_grid = 0.0:0.00000001:0.00000002;
feasibles = zeros(size(h_grid));

for i = 2:length(h_grid)
    h = h_grid(i);
    tau = 0.5 * h;
    [P, Q, tmin, xfeas, feasible] = verify_stability_1(h)
    if feasible
        break
    end
    feasibles(i) = feasible;
end

figure;
hold on;
plot(h_grid, feasibles);
grid;
ylabel("rho")
xlabel("Sampling time")


function [F, G] = get_disc_system(h)
   ogss = create_ogss();

   F = expm(ogss.A*h);
   G = (expm(ogss.A*h) - eye(size(ogss.A))) * inv(ogss.A) * ogss.B
end


function [P, Q, tmin, xfeas, feasible] = verify_stability_1(h)
    ogss = create_ogss();
    
    poles = [-2+1i -2-1i];
    K1 = place(ogss.A, ogss.B, poles);
    
    poles = [-1 -3];
    K2 = place(ogss.A, ogss.B, poles);

    [Fh1, Gh1] = get_disc_system(h);
    [Fh2, Gh2] = get_disc_system(2*h);
    
    h1update1 = Fh1 - Gh1 * K1;
    h1update2 = Fh1 - Gh1 * K2;

    h2update1 = Fh2;
    h2update2 = Fh2 - Gh2 * K2;

    evo1 = h2update1 * h1update1;
    evo2 = h2update2 * h1update1;
    evo3 = h2update1 * h1update2;
    evo4 = h2update2 * h1update2;

    setlmis([]);
    
    P = lmivar(1, [2, 1]);
    Q = lmivar(1, [2, 1]);
    
    lmiterm([1, 1, 1, P], evo1', evo1);
    lmiterm([1, 1, 1, P], -eye(2), eye(2));
    lmiterm([-1, 1, 1, Q], -eye(2), eye(2));
    
    lmiterm([2, 1, 1, P], evo2', evo2);
    lmiterm([2, 1, 1, P], -eye(2), eye(2));
    lmiterm([-2, 1, 1, Q], -eye(2), eye(2));

    lmiterm([3, 1, 1, P], evo3', evo3);
    lmiterm([3, 1, 1, P], -eye(2), eye(2));
    lmiterm([-3, 1, 1, Q], -eye(2), eye(2));
    
    lmiterm([4, 1, 1, P], evo4', evo4);
    lmiterm([4, 1, 1, P], -eye(2), eye(2));
    lmiterm([-4, 1, 1, Q], -eye(2), eye(2));
    
    lmiterm([-5, 1, 1, P], eye(2), eye(2));
    lmiterm([-6, 1, 1, Q], eye(2), eye(2));

    % 
    % newlmi();
    % lmiterm([2, 1, 1, P], evo1', evo1);
    % lmiterm([2, 1, 1, P], -eye(2), eye(2));

    lmis = getlmis();
    [tmin, xfeas] = feasp(lmis);
    P = dec2mat(lmis, xfeas, P);
    Q = 1;
    Q = dec2mat(lmis, xfeas, Q);

    feasible = sign(-tmin);
end