 %% System models
s = tf('s');
z = tf('z');
h = 0.004;  % sampling time
% plant_tf = 1 / (s * (s + 2) * (s*s + 100*s + 2600));

% controller_1 = 2.94e5*(s+2)/(s+20);
% controller_2 = 1.1e6*(s+2)*(s+2)/(s+20)/s;

plant = 1 / (s * (s + 2) * (s*s + 100*s + 2600));
plant.u = 'u';
plant.y = 'y';

controller_1 = 2.94e5*(s+2)/(s+20);
controller_1.y = 'uC';
controller_1.u = 'e';

controller_2 = 1.1e6*(s+2)*(s+2)/(s+20)/s;
controller_2.y = 'uC';
controller_2.u = 'e';

dist = sumblk('u = du + uC');
err = sumblk('e = r - y');

complete_ss = connect(controller_1, plant, err, dist, {'r', 'du'}, 'y');

%% Sampling period
sampling_frequency =  2 * max(abs(eigs(complete_ss.A)));
sampling_period = round(2 * pi / sampling_frequency, 3, 'significant', TieBreaker='tozero');
sampling_period = 2*pi/(20*7.14);
%% Discrete time system

% plant_disc = tf(plant, sampling_period);
% controller_1_disc = tf(controller_1, sampling_period);
% controller_2_disc = tf(controller_2_disc, sampling_period);
% 
% complete_ss_disc = ss(complete_ss, sampling_period);
% step(complete_ss_disc(1))

plant_disc = make_discrete(plant, sampling_period);

% controller_1_disc = make_discrete(controller_1, sampling_period);
% controller_2_disc = make_discrete(controller_2, sampling_period);

controller_1_disc = c2d(controller_1, sampling_period, 'tustin');
controller_2_disc = c2d(controller_2, sampling_period, 'tustin');

complete_disc_1 = connect(controller_1_disc, plant_disc, err, dist, {'r'}, 'y');
complete_disc_2 = connect(controller_2_disc, plant_disc, err, dist, {'du'}, 'y');

figure
step(complete_disc_1)
stepinfo(complete_disc_1)
pole(complete_disc_1)

figure
step(complete_disc_2)
stepinfo(complete_disc_2)
pole(complete_disc_2)

% tf(plant_disc)  %just checking
% plant_disc = c2d(plant, sampling_period, 'zoh')
%% Full-state-feedback controller

controllable = rank([plant_disc.B ...
    plant_disc.A*plant_disc.B ...
    plant_disc.A*plant_disc.A*plant_disc.B ...
    plant_disc.A*plant_disc.A*plant_disc.A*plant_disc.B]')

% pole set from continuous
p1_cont = [-53.15+14.91*1i -53.15-14.91i -6.85+7.04*1i -6.85-7.04*1i];
% poles that were just made up
p2_cont = [-52 -54 -14+14.2*1i -14-14.2*1i];
p3_cont = [-60+10 -60-10 -7.5+7.4*1i -7.5-7.4*1i];


[K1, tracking1] = pole_placement(plant_disc, p1_cont, sampling_period);
[fsf_y1, fsf_t1] = step(tracking1(1));
stepinfo(fsf_y1, fsf_t1)

[K2, tracking2] = pole_placement(plant_disc, p2_cont, sampling_period);
[fsf_y2, fsf_t2] = step(tracking2(1));
stepinfo(fsf_y2, fsf_t2)

[K3, tracking3] = pole_placement(plant_disc, p3_cont, sampling_period);
[fsf_y3, fsf_t3] = step(tracking3(1));
stepinfo(fsf_y3, fsf_t3)

[K4, tracking4] = pole_placement(plant_disc, p2_cont, 0.5*sampling_period);
[fsf_y4, fsf_t4] = step(tracking4(1));
stepinfo(fsf_y4, fsf_t4)

figure;
hold on

[fsf_y1, fsf_t1] = zoh(fsf_y1, fsf_t1);
[fsf_y2, fsf_t2] = zoh(fsf_y2, fsf_t2);
[fsf_y3, fsf_t3] = zoh(fsf_y3, fsf_t3);
[fsf_y4, fsf_t4] = zoh(fsf_y4, fsf_t4);

plot(fsf_t1, fsf_y1, "DisplayName", "First set of poles");
plot(fsf_t2, fsf_y2, "DisplayName", "Secondcl set of poles");
plot(fsf_t3, fsf_y3, "DisplayName", "Third set of poles");
plot(fsf_t4, fsf_y4, "DisplayName", "Double sample rate");

legend
hold off
%% Observer design
% p1_obs_cont = [-260 -270 -75+76i -75-76i];
p1_obs_cont = 5*p2_cont;

observable = rank([plant_disc.C; ...
    plant_disc.C*plant_disc.A; ...
    plant_disc.C*plant_disc.A*plant_disc.A; ...
    plant_disc.C*plant_disc.A*plant_disc.A*plant_disc.A]')

[K_obs, L_obs, observed_tracking] = observer_placement(plant_disc, p2_cont, p1_obs_cont, sampling_period);
[K_obs_drate, L_obs_drate, observed_tracking_drate] = observer_placement(plant_disc, p2_cont, p1_obs_cont, 0.5*sampling_period);


inital_state = [1 10 100 1000];
x0 = [inital_state, 0*inital_state];
observer_config = RespConfig("InitialState", x0);
tracking_config = RespConfig("InitialState", inital_state);

% step responses
[obs_y, obs_t] = step(observed_tracking, observer_config);
observer_stepinfo = stepinfo(obs_y(:, 1), obs_t)

[trk_y, trk_t] = step(tracking2(1), tracking_config);
tracker_stepinfo = stepinfo(trk_y, trk_t)

[obs_y2, obs_t2] = step(observed_tracking_drate, observer_config);
observer_stepinfo_drate = stepinfo(obs_y2, obs_t2)

[trk_y2, trk_t2] = step(tracking4(1), tracking_config);
tracker_stepinfo_drate = stepinfo(trk_y2(:, 1), trk_t2)

% plots
obs_err = obs_y(:, 5:8) - obs_y(:, 1:4);
obs_err2 = obs_y2(:, 5:8) - obs_y2(:, 1:4);
[p1y, p1t] = zoh(obs_y(:, 1), obs_t);
[p1y2, p1t2] = zoh(trk_y, trk_t);

figure;
hold on
plot(p1t, p1y, 'DisplayName', 'Output feedback')
plot(p1t2, p1y2, 'DisplayName', "State feedback")
% this is where i found out about stairs :(
stairs(obs_t2, obs_y2(:, 1), 'DisplayName', "Double rate output feedback")
stairs(trk_t2, trk_y2, 'DisplayName', "Double rate state feedback")
legend
hold off

% error of the state estimator plot
[p2e1, p2t] = zoh(obs_err(:, 1), obs_t);
[p2e2, p2t] = zoh(obs_err(:, 2), obs_t);
[p2e3, p2t] = zoh(obs_err(:, 3), obs_t);
[p2e4, p2t] = zoh(obs_err(:, 4), obs_t);

figure;
subplot(2, 2, 1)
hold on
plot(p2t, p2e1);
stairs(obs_t2, obs_err2(:, 1));
hold off
subplot(2, 2, 2)
hold on
plot(p2t, p2e2);
stairs(obs_t2, obs_err2(:, 2));
hold off
subplot(2, 2, 3)
hold on
plot(p2t, p2e3);
stairs(obs_t2, obs_err2(:, 3));
hold off
subplot(2, 2, 4)
hold on
plot(p2t, p2e4);
stairs(obs_t2, obs_err2(:, 4));
hold off

%% Revised observer design
rev_obs = p2_cont;
rev_obs_cont = 5 * rev_obs;

% unintegrated
L_int_1 = 0.00;
[K_obs, L_obs, observed_tracking] = observer_placement_integrator(plant_disc, rev_obs, rev_obs_cont, L_int_1, 0.5*sampling_period);

inital_state_int = [0 0 0 0];
x0_int = [inital_state_int, 0*inital_state_int, 0, 1];
observer_config_int = RespConfig("InitialState", x0_int);

figure;
[y_obs_int, t_obs_int] = lsim(observed_tracking(1), zeros(100, 1)', [], x0_int);
stairs(t_obs_int, y_obs_int)



% integrated
L_int_1 = 0.02;
[K_obs, L_obs, observed_tracking] = observer_placement_integrator(plant_disc, rev_obs, rev_obs_cont, L_int_1, 0.5*sampling_period);

inital_state_int = [1 10 10 1000];
x0_int = [inital_state_int, 0*inital_state_int, 0, 1];
observer_config_int = RespConfig("InitialState", x0_int);


pole(observed_tracking(1))
huh = ss(observed_tracking.A(1:9, 1:9), observed_tracking.B(1:9), observed_tracking.C(1, 1:9), observed_tracking.D(1));
pole(huh)

huh = ss(observed_tracking.A(1:8, 1:8), observed_tracking.B(1:8), observed_tracking.C(1, 1:8), observed_tracking.D(1));
pole(huh)

figure;
[y_obs_int, t_obs_int] = lsim(observed_tracking(1), zeros(100, 1)', [], x0_int);
stairs(t_obs_int, y_obs_int)

figure;
[y_obs_int, t_obs_int] = lsim(observed_tracking(1), ones(200, 1)', [], x0_int);
stairs(t_obs_int, y_obs_int)

%% Discrete lqr control
disp("LQR stepinfos")

% first set of weights
md_1 = (0.05*[1 102 2600 5200]).^2;
mu_1 = 1000^2;
Q1 = diag(1./md_1);
R1 = 1/mu_1;
N_1 = [0, 0, 0, 0]';
[K_lqr, S_lqr, P_lqr, lqr_tracking1] = dlqr_creation(plant_disc, Q1, R1, N_1, sampling_period);
[lqr_y1, lqr_t1] = step(lqr_tracking1(1));
stepinfo(lqr_y1, lqr_t1, "SettlingTimeThreshold", 0.01)

% second set of weights
md_2 = (0.05*[0.7 Inf Inf Inf]).^2;
mu_2 = 1000^2;
Q2 = diag(1./md_2);
R2 = 1/mu_2;
N_2 = [0, 0, 0, 0]';
[K_lqr2, S_lqr2, P_lqr2, lqr_tracking2] = dlqr_creation(plant_disc, Q2, R2, N_2, sampling_period);
[lqr_y2, lqr_t2] = step(lqr_tracking2(1));
stepinfo(lqr_y2, lqr_t2, "SettlingTimeThreshold", 0.01)

% third set of weights
md_3 = (0.05*[Inf Inf Inf 5200*0.7]).^2;
mu_3 = 100000^2;
Q3 = diag(1./md_3);
R3 = 1/mu_3;
N_3 = [0, 0, 0, 0]';
[K_lqr3, S_lqr3, P_lqr3, lqr_tracking3] = dlqr_creation(plant_disc, Q3, R3, N_3, sampling_period);
[lqr_y3, lqr_t3] = step(lqr_tracking3(1));
stepinfo(lqr_y3, lqr_t3, "SettlingTimeThreshold", 0.01)

% fourth set of weights
md_4 = (0.05*[1 102 2600 5200]).^2;
mu_4 = 100000^2;
Q4 = diag(1./md_4);
R4 = 1/mu_4;
N_4 = [0, 0, 0, 0]';
[K_lqr4, S_lqr4, P_lqr4, lqr_tracking4] = dlqr_creation(plant_disc, Q4, R4, N_4, sampling_period);
[lqr_y4, lqr_t4] = step(lqr_tracking4(1));
stepinfo(lqr_y4, lqr_t4, "SettlingTimeThreshold", 0.01)

% plots
figure;
hold on
stairs(lqr_t1, lqr_y1, 'DisplayName', 'First weights');
stairs(lqr_t2, lqr_y2, 'DisplayName', 'Second weights');
stairs(lqr_t3, lqr_y3, 'DisplayName', 'Third weights');
stairs(lqr_t4, lqr_y4, 'DisplayName', 'Fourth weights');
legend
hold off
%% Limited control inputs

% lqr input sizesclc
[lqr_yu, lqr_tu] = step(lqr_tracking4(2));

% Output feedback input sizes
[obs_yu, obs_tu] = step(observed_tracking_drate(9));

% Full state feedback input sizes
[fsf_yu, fsf_tu] = step(tracking2(2));

% PID input sizes
PD_usize = connect(controller_1_disc, plant_disc, err, dist, 'r', {'y', 'u'});
[pid_yu, pid_tu] = step(PD_usize(2));

[lqr_yu, lqr_tu] = zoh(lqr_yu, lqr_tu);
[obs_yu, obs_tu] = zoh(obs_yu, obs_tu);
[fsf_yu, fsf_tu] = zoh(fsf_yu, fsf_tu);
[pid_yu, pid_tu] = zoh(pid_yu, pid_tu);

figure;
hold on
plot(lqr_tu, lqr_yu, "DisplayName", "LQR input size");
plot(obs_tu, obs_yu, "DisplayName", "Output feedback input size");
plot(fsf_tu, fsf_yu, "DisplayName", "Full state feedback input size");
plot(pid_tu, pid_yu, "DisplayName", "PID input size");
legend
hold off

%% Retune PID
% controlsystemdesigner:
% reduce pole until umax is below 1000 -> 
% increase gain until damping is 0.685 -> repeat
% new_disc_zero = 0.18;
% new_disc_pole = 0.97;
% new_disc_gain = 49.391;

new_disc_zero = 0.89;
new_disc_pole = 0.95321;
new_disc_gain = 536.69;

new_pd_disc = zpk(new_disc_zero, new_disc_pole, new_disc_gain, sampling_period);
new_pd_disc.u = 'e';
new_pd_disc.y = 'u';
PD_new = connect(new_pd_disc, plant_disc, err, 'r', {'y', 'u'});
stepinfo(PD_new(1), "SettlingTimeThreshold", 0.01)
stepinfo(PD_new(2))
rlocus(PD_new(1))
% SettlingTime: 13.3320
% Peak: 994.2225
%% Retune fsf
pole_centers = -0.83857;
pole_centers = [pole_centers-0.01 pole_centers+0.01];
retune_poles = [pole_centers(1)+0.003i pole_centers(1)-0.003i 
    pole_centers(2)+0.01i pole_centers(2)-0.01i];
% SettlingTime: 12.1440

% pole_centers = -0.5104;
% pole_centers = [pole_centers-10 pole_centers+0.01];
% retune_poles = [pole_centers(1)+0.003i pole_centers(1)-0.003i 
%     pole_centers(2)+0.01i pole_centers(2)-0.01i];
% 13.4640
[K_new, tracking_new] = pole_placement(plant_disc, retune_poles, sampling_period);
stepinfo(tracking_new(1), "SettlingTimeThreshold", 0.01)
stepinfo(tracking_new(2))
step(tracking_new)
%% Retune lqr

% fourth set of poles
% md_new = [0.05 102 2600 5200];
% mu_new = 50;
% N_new = [-0.5, 0, 0, 0]';

mu = [1 102 2600 5200];
mr = 492;
Q_new1 = diag(1./(mu.*mu));
R_new1 = 1/(mr*mr);
N_new1 = [0; 0; 0; 0];
[K_lqrnew1, S_lqrnew1, P_lqrnew1, lqr_trackingnew1] = dlqr_creation(plant_disc, Q_new1, R_new1, N_new1, sampling_period);

stepinfo(lqr_trackingnew1(1), "SettlingTimeThreshold", 0.01)
stepinfo(lqr_trackingnew1(2))
step(lqr_trackingnew1)

mu = [Inf Inf Inf Inf];
mr = 10;
d1 = [1 0 1 0];
Q_new2 = [d1(1)*d1(1) d1(1)*d1(2) d1(1)*d1(3) d1(1)*d1(4);
        d1(1)*d1(2) d1(2)*d1(2) d1(2)*d1(3) 0;
        d1(1)*d1(3) d1(2)*d1(3) d1(3)*d1(3) 0;
        d1(1)*d1(4) 0 0 d1(4)*d1(4)];
Q_new2 = Q_new2 + diag(1./(mu.*mu));
R_new2 = mr*mr;

dr = [0 0 100 1 280];
N_new2 = (dr(1:4)*dr(5))';
R_new2 = R_new2 + dr(5)*dr(5);
Q_new2 = Q_new2 + diag(dr(1:4).*dr(1:4))
[K_lqrnew2, S_lqrnew2, P_lqrnew2, lqr_trackingnew2] = dlqr_creation(plant_disc, Q_new2, R_new2, N_new2, sampling_period);

stepinfo(lqr_trackingnew2(1), "SettlingTimeThreshold", 0.01)
stepinfo(lqr_trackingnew2(2))
step(lqr_trackingnew2)

%% combined plots

[pd_tx, pd_yx] = step(PD_new(1));
[pd_tu, pd_yu] = step(PD_new(2));

[trk_tx, trk_yx] = step(tracking_new(1));
[trk_tu, trk_yu] = step(tracking_new(2));

[lqr_tx, lqr_yx] = step(lqr_trackingnew2(1));
[lqr_tu, lqr_yu] = step(lqr_trackingnew2(2));

figure
hold on
stairs(pd_yx, pd_tx, "DisplayName", "PD");
stairs(trk_yx, trk_tx, "DisplayName", "Full state feedback");
stairs(lqr_yx, lqr_tx, "DisplayName", "LQR");
legend
hold off

figure

hold on
stairs(pd_yu, pd_tu, "DisplayName", "PD");
stairs(trk_yu, trk_tu, "DisplayName", "Full state feedback");
stairs(lqr_yu, lqr_tu, "DisplayName", "LQR");
legend
hold off

%% PD disturbance
new_disc_zero = 0.89;
new_disc_pole = 0.95321;
new_disc_gain = 477.36;

Li = 0.004;

new_pd_disc_int = new_pd_disc * (1+Li/(z-1));
new_pd_disc_int.u = 'e';
new_pd_disc_int.y = 'uC';

PD_new = connect(new_pd_disc, plant_disc, err, 'r', {'y', 'u'});
stepinfo(PD_new(1), "SettlingTimeThreshold", 0.01)
stepinfo(PD_new(2))
rlocus(PD_new(1))

PD_new_int = connect(new_pd_disc_int, plant_disc, err, dist, {'r', 'du'}, {'y', 'u'});
stepinfo(PD_new_int(1, 1), "SettlingTimeThreshold", 0.01)
stepinfo(PD_new_int(1, 2))

[pid_dist_yy, pid_dist_ty] = step(PD_new_int(1, 1));
[pid_dist_yu, pid_dist_tu] = step(PD_new_int(2, 1));
[pid_dist_yd, pid_dist_td] = step(PD_new_int(1, 2));


figure;
stairs(pid_dist_td, pid_dist_yd)
figure;
stairs(pid_dist_tu, pid_dist_yu)
figure;
stairs(pid_dist_ty, pid_dist_yy)

%% Disturbance observer

[K_disturbance_observer, L_disturbance_observer, disturbance_observer] = observer_placement_disturbance(plant_disc, p2_cont, p1_obs_cont, 0.5*sampling_period);
inital_state_dist = [0 0 0 0 1];
x0 = [inital_state_dist, 0*inital_state_dist];
observer_config = RespConfig("InitialState", x0);
tracking_config = RespConfig("InitialState", inital_state_dist);

% step responses
% [obs_dy, obs_dt] = step(disturbance_observer(1), zeros());
% observer_stepinfo = stepinfo(obs_y(:, 1), obs_t);

% lsim(disturbance_observer(1), zeros(100, 1)', [], x0);
[obs_dy, obs_dt] = lsim(disturbance_observer([1, 11]), zeros(100, 1)', [], x0);

figure;
stairs(obs_dt, obs_dy(:, 1));
figure; 
stairs(obs_dt, obs_dy(:, 2));


%% State integrator
pole_centers = -0.83857;
pole_centers = [pole_centers-0.01 pole_centers+0.01];
retune_poles = [pole_centers(1)+0.003i pole_centers(1)-0.003i 
    pole_centers(2)+0.01i pole_centers(2)-0.01i];

pole_center = -0.4;
distance = 0.03;
arm1 = distance+distance*1i;
dist_poles = [pole_center, pole_center+arm1, pole_center+arm1*1i, pole_center-arm1, pole_center-arm1*1i]';
[K_state_integrator, integ_tracking] = pole_placement_disturbance(plant_disc, dist_poles, sampling_period);

x0_integ = [0 0 0 0 0 1];
% lsim(integ_tracking, zeros(1000, 1)', [], x0_integ);

figure;
% step(integ_tracking(1))

% [obs_dy, obs_dt] = lsim(integ_tracking([1, 2]), zeros(500, 1)', [], x0);
[fsf_yd, fsf_td] = lsim(integ_tracking(1), zeros(1, 1000), [], x0_integ);
[fsf_yu, fsf_tu] = step(integ_tracking(2));
[fsf_yy, fsf_ty] = step(integ_tracking(1));

figure;
stairs(fsf_td, fsf_yd)
figure;
stairs(fsf_tu, fsf_yu)
figure;
stairs(fsf_ty, fsf_yy)


%% lqr integrator
[K_lqrint, S_lqrint, P_lqrint, lqr_trackingint] = dlqr_creation_disturbance(plant_disc, Q_new2, R_new2, N_new2, sampling_period);


x0_integ = [0 0 0 0 0 1];
% lsim(integ_tracking, zeros(1000, 1)', [], x0_integ);

figure;
% step(integ_tracking(1))

% [obs_dy, obs_dt] = lsim(integ_tracking([1, 2]), zeros(500, 1)', [], x0);
[lqr_yd, lqr_td] = lsim(lqr_trackingint(1), zeros(1, 1000), [], x0_integ);
[lqr_yu, lqr_tu] = step(lqr_trackingint(2));
[lqr_yy, lqr_ty] = step(lqr_trackingint(1));

figure;
stairs(lqr_td, lqr_yd)
figure;
stairs(lqr_tu, lqr_yu)
figure;
stairs(lqr_ty, lqr_yy)


%% Delay in the system



%% Functions

function [sys_disc] = make_discrete(sys, h)
    if isa(sys, "tf")
        sys = compreal(minreal(sys), 'c');
        % canonical form
        inverter =  fliplr(eye(size(sys.A)));
        sys = ss2ss(sys, inverter);
    end
    bk = size(sys.B);
    n = bk(1);
    % sys.A
    % det(sys.A)
    % assert(det(sys.A) ~= 0);
    phi = expm(sys.A*h);
    psi = integral(@(t) expm(sys.A.*t), 0, h, 'ArrayValued', true);
    gamma = psi * sys.B;

    sys_disc = ss(phi, gamma, sys.C, sys.D, h);
    sys_disc.u = sys.u;
    sys_disc.y = sys.y;
end

function [K, tracking] = pole_placement(plant_disc, p_cont, sampling_period)
    A = plant_disc.A;
    B = plant_disc.B;
    C = plant_disc.C;
    D = plant_disc.D;
    p = exp(p_cont*sampling_period);
    K = place(A, B, p);
    
    
    placed = ss(A-B*K, B, C, D, sampling_period);
    comb_C = [C; -K];
    comb_D = [0; 1/dcgain(placed)];
    tracking = ss(A-B*K, B/dcgain(placed), comb_C, comb_D, sampling_period);
    
    pole(placed);
    % damp(placed);

    tracking_gain = dcgain(placed);
end

function [K, L, tracking] = observer_placement(plant_disc, p_cont, p_cont_obs, sampling_period)
    A = plant_disc.A;
    B = plant_disc.B;
    C = plant_disc.C;
    D = plant_disc.D;
    disp("new observer placement: ")
    p = exp(p_cont*sampling_period);
    p_obs = exp(p_cont_obs*sampling_period);
    
    K = place(A', C', p_obs)';
    L = place(A, B, p);

    comb_A = [A, -B*L; K*C, A-K*C-B*L];
    comb_B = [B; B];
    comb_C = [C zeros(1, 4)];
    comb_D = 0;

    observed = ss(comb_A, comb_B, comb_C, comb_D, sampling_period);
    
    output_C = eye(9, 8);
    output_C(9, 5:8) = -L;
    tracking = ss(comb_A, comb_B/dcgain(observed), output_C, flip(eye(9, 1))/dcgain(observed), sampling_period);

    pole(observed)
    damp(observed)
    if any(imag(K) ~= 0)
        disp("!K has imaginary values!")
        K
    end
    if any(imag(L) ~= 0)
        disp("!L has imaginary values!")
        L
    end
    
    observed_gain = dcgain(observed);
end

function [K, S, P, tracking] = dlqr_creation(plant_disc, Q, R, N, sampling_period)
    A = plant_disc.A;
    B = plant_disc.B;
    C = plant_disc.C;
    D = plant_disc.D;
    
    [K, S, P] = dlqr(A, B, Q, R, N);
    lqed = ss(A-B*K, B, C, D, sampling_period);
    comb_C = [C; -K];
    comb_D = [0; 1/dcgain(lqed)];
    tracking = ss(A-B*K, B/dcgain(lqed), comb_C, comb_D, sampling_period);

    lqr_gain = dcgain(lqed);
end

function [y, t] = zoh(y, t)
    y = repelem(y, 2);
    y = y(1:end-1);
    t =  repelem(t, 2);
    t = t(2:end);
end

%% disturbance models
function [K, tracking] = pole_placement_disturbance(plant_disc, p_cont, sampling_period)
    A = plant_disc.A;
    B = plant_disc.B;
    C = plant_disc.C;
    D = plant_disc.D;
    p = exp(p_cont*sampling_period);
    

    A_big = [A zeros(4, 1); C 1];
    B_big = [B; 0];

    K = place(A_big, B_big, p);
    L_int = K(5);
    K = K(1:4);
    placed = ss(A-B*K, B, C, D, sampling_period);


    L_int = L_int;
    comb_A = [A-B*K, -B*L_int, B; C, 1, 0; zeros(1, 5) 1];
    comb_B = [B; -dcgain(placed); 0];
    comb_C = [C zeros(1, 2); -K -L_int, 0; 0 0 0 0 1 0];
    comb_D = [0; 1/dcgain(placed); 0];

    tracking = ss(comb_A, comb_B/dcgain(placed), comb_C, comb_D, sampling_period);
    
    pole(placed);
    % damp(placed);

    tracking_gain = dcgain(placed);
end

function [K, L, tracking] = observer_placement_disturbance(plant_disc, p_cont, p_cont_obs, sampling_period)
    A = plant_disc.A;
    B = plant_disc.B;
    C = plant_disc.C;
    D = plant_disc.D;
    disp("new observer placement: ")

    p = exp(p_cont*sampling_period);
    p_cont_obs = [p_cont_obs -10];
    p_obs = exp(p_cont_obs*sampling_period);
    
    A_big = [A B; 0 0 0 0 1];
    B_big = [B; 0];
    C_big = [C 0];
    D_big = D;

    K = place(A_big', C_big', p_obs)';
    L = place(A, B, p);
    L = [L 1];

    comb_A = [A_big, -B_big*L; K*C_big, A_big-K*C_big-B_big*L];
    comb_B = [B_big; B_big];
    comb_C = [C_big zeros(1, 5)];
    comb_D = 0;

    observed = ss(comb_A, comb_B, comb_C, comb_D, sampling_period);
    
    output_C = eye(11, 10);
    output_C(11, 6:10) = -L;
    tracking = ss(comb_A, comb_B/dcgain(observed), output_C, flip(eye(11, 1))/dcgain(observed), sampling_period);

    pole(observed)
    damp(observed)
    if any(imag(K) ~= 0)
        disp("!K has imaginary values!")
        K
    end
    if any(imag(L) ~= 0)
        disp("!L has imaginary values!")
        L
    end
    
    observed_gain = dcgain(observed);
end

function [K, S, P, tracking] = dlqr_creation_disturbance(plant_disc, Q, R, N, sampling_period)
    A = plant_disc.A;
    B = plant_disc.B;
    C = plant_disc.C;
    D = plant_disc.D;
    
    [K, S, P] = dlqr(A, B, Q, R, N);
    lqed = ss(A-B*K, B, C, D, sampling_period);
    L_int = 0.004/dcgain(lqed);

    comb_A = [A-B*K, -B*L_int, B; C, 1, 0; zeros(1, 5) 1];
    comb_B = [B; -dcgain(lqed); 0];
    comb_C = [C zeros(1, 2); -K -L_int, 0; 0 0 0 0 1 0];
    comb_D = [0; 1/dcgain(lqed); 0];

    tracking = ss(comb_A, comb_B/dcgain(lqed), comb_C, comb_D, sampling_period);

    lqr_gain = dcgain(lqed);
end

function [K, L, tracking] = observer_placement_integrator(plant_disc, p_cont, p_cont_obs, L_int, sampling_period)
    A = plant_disc.A;
    B = plant_disc.B;
    C = plant_disc.C;
    D = plant_disc.D;
    disp("new observer placement: ")

    p = exp(p_cont*sampling_period);
    p_obs = exp(p_cont_obs*sampling_period);

    K = place(A', C', p_obs)';
    L = place(A, B, p);
    
    comb_A = [A, -B*L; K*C, A-K*C-B*L];
    comb_B = [B; B];
    comb_C = [C zeros(1, 4)];
    comb_D = 0;

    undisturbed = ss(comb_A, comb_B, comb_C, comb_D, sampling_period);
    
    L_int = L_int/dcgain(undisturbed);

    comb_A = [A, -B*L,  -B*L_int, B; 
        K*C, A-K*C-B*L, -B*L_int, zeros(4, 1);
        C zeros(1, 4) 1 0;
        zeros(1, 9) 1];
    comb_B = [B; B; -1; 0];
    comb_C = [C zeros(1, 6)];
    comb_D = 0;

    observed = ss(comb_A, comb_B, comb_C, comb_D, sampling_period);
    
    output_C = eye(10);
    tracking = ss(comb_A, comb_B/dcgain(observed), output_C, D, sampling_period);

    pole(observed)
    damp(observed)

    observed_gain = dcgain(observed);
end
