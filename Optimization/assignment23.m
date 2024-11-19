%parameters
pars = dictionary;
pars("N") = 5; %number of sections
pars("nt") = 120; %number of timesteps
pars("E1") = 15;
pars("E2") = 7;
pars("E3") = 8;
pars("rhom") = 120;                  % Max density [veh/km/lane]
pars("rhoc") = 33.5 + pars("E1")/3;  % Critical density [veh/km/lane]
pars("T") = 10/3600;                 % Time step [s]
pars("L") = 1;                       % Segment length [km]
pars("lambda") = 3;                  % Number of lanes
pars("tau") = 10/3600;               % Relaxation time [s]
pars("mu") = 80;                     % Convection term coefficient [km^2/h]
pars("Cr") = 2000;                   % Ramp capacity [veh/h]
pars("alpha") = 0.1;                 % Driver compliance factor
pars("K") = 10;                      % Model parameter [veh/km/lane]
pars("a") = 2;                       % Model parameter
pars("vf") = 110;                    % Free flow speed [km/h]
pars("Dr") = 1500;                   % Ramp demand [veh/h]
q0 = [ repmat(7000 + 100*pars("E2"), 60, 1); repmat(2000 + 100*pars("E3"), 60, 1) ];

% precompute some numbers
pars("ToverlambdaL") = pars("T") / pars("lambda") / pars("L");
pars("Tovertau") = pars("T") / pars("tau");
pars("ToverL") = pars("T") / pars("L");
pars("muTovertauL") = pars("mu") * pars("T") / pars("tau") / pars("L");

% Initial Conditions
rho_init = 30;  % Initial density for all segments [veh/km/lane]
v_init = 80;  % Initial speed for all segments [km/h]
wr_init = 0;  % Initial ramp queue length [veh]
% TTS = zeros(pars("nt"), 1);  % Total Time Spent (objective to minimize)
% rho = rho_init * ones(pars("N"), pars("nt"));
% v = v_init * ones(pars("N"),pars("nt"));
% wr = ones(pars("nt"), 1) * wr_init;
% r = 0;


%% NO CONTROL

x0 = [rho_init, rho_init, rho_init, rho_init, rho_init, ...
    v_init, v_init, v_init, v_init, v_init, wr_init];
u = repmat([120, 1], 120, 1);

x = simulate(pars, q0, u, x0);

[TTS_t, TTS_tot] = TTS_fnc(x, pars);

%% CONTROL

nonlcon = @(x) cons_fnc(x, pars, q0, x0);
[lb, ub] = create_bounds(pars);

initial_guess = repmat([60, 0.1], 120, 1);
objective = @(x) the_goal(x, pars, q0, x0);

options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp','OptimalityTolerance', 1e-2);
% [X_opt, TTS_opt, ~, ~, ~, ~, lambda_out] = fmincon(objective, initial_guess, [],[], [],[],[], [], nonlcon);
u_results = fmincon(objective, initial_guess, [],[], [],[], lb ,ub , nonlcon, options);
%% PLOTTING

x_results = simulate(pars, q0, u_results, x0);
% output conversion
rho_out = x_results(:, 1:5).';
v_out = x_results(:, 6:10).';
wr_out = x_results(:, 11).';

time = (0:pars("nt")-1) * pars("T") / 60;  % Time in minutes

% Plot density (rho) for all segments
figure;
for i = 1:pars("N")
    plot(time, rho_out(i, :));
    hold on;
end
title('Density (rho) over Time');
xlabel('Time (minutes)');
ylabel('Density (veh/km/lane)');
legend('Segment 1', 'Segment 2', 'Segment 3', 'Segment 4', 'Segment 5');
hold off;

% Plot speed (v) for all segments
figure;
for i = 1:5
    plot(time, v_out(i, :));
    hold on;
end
title('Speed (v) over Time');
xlabel('Time (minutes)');
ylabel('Speed (km/h)');
legend('Segment 1', 'Segment 2', 'Segment 3', 'Segment 4', 'Segment 5');
hold off;

% Plot ramp queue length (wr) over time
figure;
plot(time, wr_out);
title('Ramp Queue Length (wr) over Time');
xlabel('Time (minutes)');
ylabel('Queue Length (veh)');

%% STATE DYNAMICS
function [x_new] = state_update(x_old, pars, u, q0, k)
    % just variables
    x_new = zeros(11, 1);
    % rho = [0 , x_old(1:5), x_old(5)]; 
    % v = [x_old(5), x_old(6:10), 0];
    % q = [q0(k), x_old(1:5)*x_old(6:10)*pars("lambda")];
    v = x_old(6:10);
    rho = x_old(1:5);
    q = pars("lambda")*v.*rho;
    % rho(0) = q(0)/(v(0)*pars("lambda"));
    %Flow on onramp
    c1 = u(2)*pars("Cr");
    c2 = pars("Dr") + (x_old(11) / pars("T"));
    c3 = pars("Cr") * ((pars("rhom") - rho(4)) / (pars("rhom") - pars("rhoc")));
    qr4 = min([c1 c2 c3]);
    % next ramp length
    x_new(11) = x_old(11) + pars("T")*(pars("Dr")-qr4);
 
    % density update
    x_new(1) = rho(1) + pars("ToverlambdaL")*(q0(k) - q(1));
    x_new(2) = rho(2) + pars("ToverlambdaL")*(q(1) - q(2));
    x_new(3) = rho(3) + pars("ToverlambdaL")*(q(2) - q(3));
    x_new(4) = rho(4) + pars("ToverlambdaL")*(q(3) - q(4) + qr4);
    x_new(5) = rho(5) + pars("ToverlambdaL")*(q(4) - q(5));
    
    % desired speed
    bigv = [ min([ (1+pars("alpha")) * 120, ...
        pars("vf") * exp( - (rho(1) / pars("rhoc"))^pars("a") / pars("a"))]);
        min([ (1+pars("alpha")) * u(1), ...
        pars("vf") * exp( - (rho(2) / pars("rhoc"))^pars("a") / pars("a"))]);
        min([ (1+pars("alpha")) * u(1), ...
        pars("vf") * exp( - (rho(3) / pars("rhoc"))^pars("a") / pars("a"))]);
        min([ (1+pars("alpha")) * 120, ...
        pars("vf") * exp( - (rho(4) / pars("rhoc"))^pars("a") / pars("a"))]);
        min([ (1+pars("alpha")) * 120, ...
        pars("vf") * exp( - (rho(5) / pars("rhoc"))^pars("a") / pars("a"))]); ];
    
    %anticipation term
    ant = [ (pars("muTovertauL") * (rho(2) - rho(1))) / ((rho(1) + pars("K")))
        (pars("muTovertauL") * (rho(3) - rho(2))) / ((rho(2) + pars("K")))
        (pars("muTovertauL") * (rho(4) - rho(3))) / ((rho(3) + pars("K")))
        (pars("muTovertauL") * (rho(5) - rho(4))) / ((rho(4) + pars("K")))
        0];

    % next speed
    x_new(6) = v(1) + ( pars("Tovertau") * (bigv(1) - v(1)) ) + pars("ToverL") * v(1) * (v(1) - v(1)) - ant(1);
    x_new(7) = v(2) + ( pars("Tovertau") * (bigv(2) - v(2)) ) + pars("ToverL") * v(2) * (v(1) - v(2)) - ant(2);
    x_new(8) = v(3) + ( pars("Tovertau") * (bigv(3) - v(3)) ) + pars("ToverL") * v(3) * (v(2) - v(3)) - ant(3);
    x_new(9) = v(4) + ( pars("Tovertau") * (bigv(4) - v(4)) ) + pars("ToverL") * v(4) * (v(3) - v(4)) - ant(4);
    x_new(10) = v(5) + ( pars("Tovertau") * (bigv(5) - v(5)) ) + pars("ToverL") * v(5) * (v(4) - v(5)) - ant(5);
end

% this function signature is stupid at in only simulates for 119 timesteps
function [x] = simulate(pars, q0, u, x0)
    x = zeros(pars("nt"), 13);
    x(1, 1:11) = x0;
    x(1, 12:13) = u(1, :);
    for i = 2:pars("nt")
        x(i, 1:11) = state_update(x(i-1, 1:11), pars, x(i-1, 12:13), q0, i-1);
        x(i, 12:13) = u(i, :);
    end
end

function [TTS_t, TTS_tot] = TTS_fnc(x, pars)
    TTS_t = zeros(pars("nt"), 1);
    for i = 1:pars("nt")
        TTS_t(i) = pars("T") * x(i, 11) + pars("T") * pars("L") * pars("lambda") * sum(x(i, 1:5));
    end
    TTS_tot = sum(TTS_t);
end

%% FMINCON FUNCTIONS

function TTS = the_goal(u, pars, q0, initial_state)
    x = simulate(pars, q0, u, initial_state);
    [~, TTS] = TTS_fnc(x, pars);
end

function [c, ceq] = cons_fnc(u, pars, q0, initial_state)
    x = simulate(pars, q0, u, initial_state);
    c = x(:, 5:10) - pars("rhoc");
    ceq = 0;
end

% function [c, ceq] = alt_cons_fnc(u, pars, q0, initial_state)
%     c = -1;
%     x_prediction = simulate(pars, q0, u, initial_state);
%     for i = 2:pars("nt")
%         x_prediction(i, :) = state_update(x(i-1, 1:11), pars, x(i-1, 12:13), q0, i);
%     end
%     ceq = x(:, 1:11) - x_prediction;
% end

function [lb, ub] = create_bounds(pars)
    zero = zeros(pars("nt"), 1);
    one = ones(pars("nt"), 1);

    % lb
    lb_vsl = 60*one;
    lb_r = zero;
    lb = [lb_vsl; lb_r];

    % ub
    ub_vsl = 120*one;
    ub_r = one;
    ub = [ub_vsl; ub_r];
end