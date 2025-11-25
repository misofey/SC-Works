close all
% States: omega, z1_dot, z1, z2_dot, z2
% Inputs: beta, tau(generator torque), wind_speed
ops =  bodeoptions;
s =tf('s');
mimo = ss(A, B, C, D);
mimo.u = {'Beta', 'Tau', 'Wind Speed'};
mimo.y = {'omega', 'z'};

mimo = tf(mimo);
mimo = minreal(-mimo);

%%% Controller analysis %%%
bode(mimo)
omega = 0.3*2*pi;

G_0 = evalfr(mimo, omega);
G_0_inv_T = pinv(evalfr(mimo, omega))';

RGA = G_0.*G_0_inv_T

%%% Controller design %%%
% Low frequency attunation
a = 10e-4; 
% H infinity norm for sensitivity
M = 3;
% Cut off frequency
omega_b = 0.3*2*pi;

% Controller weights
W_p1 = inv(tf(makeweight(10^-4, [omega_b, 0.708], 3)));

%W_p1 =  inv(10^(-4)*(s/(2.587*10^(-4)) + 1)/(s/7.761 + 1));
W_p2 = 0.05;
W_p = [W_p1 0;
       0 W_p2];
W_p_inv = [inv(W_p1) 0;
       0 inv(W_p2)];

W_u1 = 0.005;
W_u2 = (5e-3*s^2 + 7e-4*s + 5e-5)/(s^2 + 14e-4*s + 1e-6); 
W_u = [W_u1 0    0;
       0    W_u2 0 ;
       0    0    0];
W_u_inv = [inv(W_u1) 0;
       0    inv(W_u2) ];


%% Mixed-sensitivity design %%
P = augw(mimo, W_p, W_u, []);
[K,CL,gamma] = hinfsyn(P,2,3);
minreal(K);
L = tf(minreal(mimo*K));
I = eye(size(L));
S = feedback(I,L); 
T = I-S;

bode(L)
figure
getGainCrossover(L(1,1), 1)
margin(L(1, 1))
figure
getGainCrossover(L(2,1), 1)
margin(L(2,1))
figure

ops.Title.String = "Sensitivity";
bode(S, W_p_inv)
legend(["Sensitivity", "Target Sensitivity"])
grid
figure
bode(K*S, W_u_inv)
legend(["Control Effort", "Target Control Effort"])
grid
bode(L)
grid
figure
sigma(CL, ss(GAM))
legend(["Close loop singular values", "Max singular value"])
grid
figure
nyquist(CL)
figure
%% Simulation %%
step(T(1:2, 1))
grid
figure
step(S)
grid
figure

info = stepinfo(T);
info(1, 1)
info(2, 1)

