function lumped_param_sim_wrap()

%Parameters go here
T_0 = 293;
kappa = 17e-5;
sigma_0 = 5.8e7;
delta_sigma = 0;
L_rl = 1.0005;
L_rh = 1.005;
L_0 = 1;
R_0 = 2;
L_roll = 1;
Area_rod = 0.0001;
R_L = 1;
R_H = 3;
K_fan = 0.145;
K_m = 1;
b_fan = 05;
epsilon_b = 5e9;
phi_0 = 1.55e-5;
V_fuel = 1;
rho_A = 1;
V_F = 0.027;
C_A = 1000;
R_duct = 1.44;
R_AF = 0.004;
C_F = 84410;
R_F0 = 0.01/sqrt(0.5);
V = 10;
J = 0.005;
L_motor = 0.02;

i_th_lo = 4;
i_th_hi = 4.55;

%Convert the transfer function to a state-space model
b = 1/6;
a = 1/6;
tf_rod = tf(b, [1, a]);
sys_rod = ss(tf_rod);
A = sys_rod.A;
B = sys_rod.B;
C = sys_rod.C;

x0 = [0; 0; T_0; T_0; T_0/C]; %TODO Define;

dt = 0.001;
sim_time = 2000; %[s]
steps = sim_time/dt;
t = linspace(0, sim_time, steps);
%x has all of the states over time
x = zeros(5, steps);
x(:, 1) = x0; %Initial condition


for k = 1:steps-1
    x_dot = state_derivative(x(:, k));
    if isnan(x_dot(1))
        x(1, k+1) = 0;
        x(2:5, k+1) = x(2:5, k) + x_dot(2:5)*dt;
    else
        x(:, k+1) = x(:, k) + x_dot*dt;
    end
end

%Because the state space model has a C != 1, convert to the actual value so
%that it can be plotted
T_r_allsteps = x(5, :) .* C;

%Different plot for each figures as the units are different and different
%scales
tiledlayout(3, 2)
nexttile
plot(t, x(1, :));
title("Circuit current");
xlabel("Time [s]");
ylabel("Current [A]");

nexttile
plot(t, x(2, :));
title("Motor angular velocity");
xlabel("Time [s]");
ylabel("Angular velocity [rad/s]");

nexttile
plot(t, x(3, :));
title("Furnace air temperature");
xlabel("Time [s]");
ylabel("Temeprature [K]");

nexttile
plot(t, x(4, :));
title("Furnace wall temperature");
xlabel("Time [s]");
ylabel("Temeprature [K]");


% nexttile
% plot(t, x(6, :));
% title("T_r state");
% xlabel("Time [s]");

nexttile
plot(t, T_r_allsteps);
title("Rod temeperature");
xlabel("Time [s]");
ylabel("Temperature [K]");



function x_dot = state_derivative(x)
    i = x(1);
    w = x(2);
    T_a = x(3);
    T_f = x(4);
    T_r_state = x(5);

    T_r = C*T_r_state;
    T_r_state_dot = A*T_r_state + B*(T_f);
    
    deltaT_r = T_r - T_0;
    L = L_0 * (1+kappa*deltaT_r);

    sigma = sigma_0 + delta_sigma*(T_0 - T_r)/(T_0);

    %Set value of R
    if L < L_rl
        R = inf; %This case is handled outside the loop
    elseif L < L_rh
        R = R_0 + (L-L_roll)/(sigma*Area_rod) + R_L + (R_H-R_L)*(L_rh - L)/(L_rh - L_rl);
    else
        R = R_0 + (L_rh - L_rl)/(sigma*Area_rod);
    end

    delta_P = abs(K_fan * w);
    
    %The controller based on current
    if i > i_th_hi
        V = 0;
    elseif i > i_th_lo
        V = 10;
    else
        V = 20;
    end
    
    i_dot = (V - R*i - K_m*w)/L_motor;
    w_dot = (K_m*i - b_fan * w )/J;
    T_a_dot = (epsilon_b*phi_0*V_fuel)/(rho_A*V_F*C_A) - sqrt(delta_P)*(T_a - T_0)/(R_duct*V_F) - (T_a - T_f)/(rho_A*V_F*C_A*R_AF);
    T_f_dot = (T_a - T_f)/(C_F*R_AF) - (T_f - T_0)/(C_F*R_F0);

    x_dot = [i_dot; w_dot; T_a_dot; T_f_dot; T_r_state_dot];
        
end

end

lumped_param_sim_wrap()