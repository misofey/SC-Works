% copper properties
kappa = 17e-6;  % 
gamma = 400;  %[Wm^-1K^-1]
c_r = 385;
rho = 8935;  % [kg/m^3]

%potassium properties
% kappa = 

l = 1;
area = 0.5;
N = 99;  % this N is for the points that are simulated in the heat equation, not the actual number of points in the discretization
dx = l/(N+1);

nt = 150000;
dt = 0.005;

t = nt*dt;

l_roll = 0.49;

rh = 500;
rl = 10;
lrh = 3.1;
lrl = 3.05;

sigma0 = 5;
deltasigma = 5;

t0 = 0;
t_f = 1; % furnace temperature

% temperature simulation
u = ones(N, 1)*t0;
v1 = -2 / (dx*dx) * ones(N, 1);
v2 = 1 / (dx*dx) * ones(N-1, 1);
A = diag(v1, 0);
A = A + diag(v2, 1);
A = A + diag(v2, -1);
A(N, N) = -1 / (dx*dx);
evo_matrix = eye(N) + dt*gamma/c_r/rho*A;
f = zeros(N, 1);
f(1) = t_f / dx /dx;

% simulation logging
u_hist = zeros(nt, N+2);

for i = 1:nt
    u = evo_matrix * u + dt*gamma/c_r/rho*f;
    u_hist(i, 2:N+1) = u;
    if anynan(u)
        disp("nan encountered")
    end
end
u_hist(:, 1) = t_f;
u_hist(:, N+2) = u_hist(:, N+1);
u_compl = [t_f; u; u(end)];
els = zeros(N+1, 1);

for i = 2:N
    els(i) = (u(i-1) + u(i))* kappa*dx/2 + (1-kappa*t0)*dx;
end
els(1) = (t_f + u(1))* kappa*dx/2 + (1-kappa*t0)*dx;
els(N+1) = (u(N) + u(N))* kappa*dx/2 + (1-kappa*t0)*dx;
els_compl = [0; els]; % elongations_complete
disp = cumsum(els_compl); % displacement

res_rod = 0;
i_roll = ceil((disp(end) + l - l_roll) / dx);

elf = @(temp) (1+kappa*(temp - t0));
resist = @(temp) 1/(sigma0 + deltasigma*(t0-temp/t0));
for i = i_roll:N
    r_rod = res_rod + resist(u(i-1)) * elf(u(i-1)) * dx/2 + resist(u(i)) * elf(u(i)) * dx/2;
end

if disp(end) + l < lrl
    r_res = Inf;
elseif (lrl <= (disp(end) + l)) && ((disp(end) + l) < lrh)
    r_res = (disp(end) - lrl) / (lrh - lrl) * (rh - rl) + rl;
end
% plot(u_hist(100, :))
% plot(displ)
%% curve fitting
target_value = 1-(1/exp(1));
target_line = ones(nt, 1) * target_value;
u_means = mean(u_hist(:, 2:end-2) + u_hist(:, 3:end-1), 2);
t = nt *dt;
times = linspace(0, t, nt);
lnu = -(u_means - 1);
lnln = log(lnu);
hold on;
plot(times, u_means, DisplayName="Rod temperature");
plot(times, target_line, DisplayName="Target temperature");
legend;
hold off;