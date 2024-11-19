global N E1 r E2 E3 pm rho_c T L lambda tau mu Cr rho_m alpha K a vf Dr

% Initial Conditions
rho_init = 30;  % Initial density for all segments [veh/km/lane]
v_init = 80;  % Initial speed for all segments [km/h]
wr_init = 0;  % Initial ramp queue length [veh]
time_steps = 120;
TTS = zeros(time_steps, 1);  % Total Time Spent (objective to minimize)
rho = rho_init * ones(N, time_steps);
v = v_init * ones(N,time_steps);
wr = ones(time_steps, 1) * wr_init;
% Model Parameters
E1 = 15;
E2 = 7;  % student ID parameter
E3 = 8;  % student ID parameter
pm = 120; %model paramater
rho_c = 33.5 + ((1/3) * E1);  % Critical density [veh/km/lane]
%VSL = 120;
T = 10/3600;  % Time step [s]
L = 1;  % Segment length [m]
lambda = 3;  % Number of lanes
tau = 10/3600;  % Relaxation time [s]
mu = 80;  % Convection term coefficient [km^2/h]
Cr = 2000/3600;  % Ramp capacity [veh/h]
rho_m = 120;  % Max density [veh/km/lane]
alpha = 0.1;  % Driver compliance factor
K = 10;  % Model parameter [veh/km/lane]
a = 2;  % Model parameter
vf = 110;% Free flow speed [km/h]
Dr = 1500/3600;  % Ramp demand [veh/h]
r = 0;
% Number of segments
N = 5;

% Define the optimization problem
objective = @(VSL) objective_function(N, TTS, rho, v, wr, VSL);

options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
[VSL_opt, TTS_opt, ~, ~, ~, ~, lambda_out] = fmincon(objective, 120, [],[], [],[],60,120);

disp("OPtimal VSl = " + VSL_opt)
disp("Optimal TTS = " + TTS_opt)

[~, rho_out, v_out, wr_out] = objective_function(N, TTS, rho, v, wr, VSL_opt);

time = (0:time_steps-1) * T / 60;  % Time in minutes

% Plot density (rho) for all segments
figure;
for i = 1:N
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


function [J, rho_out, v_out, wr_out] = objective_function(N, TTS, rho, v, wr, VSL)
    global  E2 E3  pm rho_c r T L lambda tau mu Cr  alpha K a vf Dr
    time_steps = 120;  % Number of simulation steps
    
    q = zeros(N,1);

    for k = 1:time_steps-1

        %Flow on onramp
        qr4 = min([r*Cr Dr+(wr(k)/T) Cr*((pm - rho(4,k))/(pm-rho_c))]);
        wr(k+1) = wr(k) + T*(Dr-qr4);
    
        for i = 1:5
            %traffic flow
            q(i) = rho(i,k)*lambda*v(i,k);
    
            if i == 1 % special condition for segement 1
                if k < 60
                    qzero = 5000 + (100*E2);
                else
                    qzero = 2000 + (100*E3);
                end
                rho(i,k+1) = rho(i,k) + ((T/(lambda*L))*(qzero - q(i))) ;
            elseif i == 4 % special condition for segement 4 including on ramp
                rho(i,k+1) = rho(i,k) + ((T/(lambda*L))*(q(i-1) - q(i)+ qr4)) ;
            else %235
                rho(i,k+1) = rho(i,k) + ((T/(lambda*L))*(q(i-1) - q(i))) ;
            end
            rho(i,k+1) = max(0, min(rho(i,k+1), rho_c));  % Keep between 0 and rho_m
            
            if i == 2 || i == 3
                Vi = min([((1+alpha)*VSL) vf*exp((-1/a)*((rho(i,k)/rho_c)^a))]);
            else
                Vi = min([((1+alpha)*120) vf*exp((-1/a)*((rho(i,k)/rho_c)^a))]);
            end
            %velocity
            if i == 5 % p6 = p5
                v(i,k+1) = v(i,k)+((T/tau)*(Vi-v(i,k)))+((T/L)*v(i,k)*(v(i-1,k)-v(i,k)));
            elseif i == 1 %v1 = v0
                v(i,k+1) = v(i,k)+((T/tau)*(Vi-v(i,k))) - ((mu*T*(rho(i+1,k)-rho(i,k)))/(tau*L*(rho(i,k)+K)));
            else
                v(i,k+1) = v(i,k)+((T/tau)*(Vi-v(i,k)))+((T/L)*v(i,k)*(v(i-1,k)-v(i,k))) - ((mu*T*(rho(i+1,k)-rho(i,k)))/(tau*L*(rho(i,k)+K)));
            end   
            v(i,k+1) = max(0, v(i,k+1));  % Velocity cannot be negative
        end  
        TTS(k) = T *wr(k) + T*sum(L * lambda * rho(:,k));
    end
    % Calculate TTS at this step
    rho_out = rho;
    v_out = v;
    wr_out = wr;
    J = sum(TTS);

    
end




