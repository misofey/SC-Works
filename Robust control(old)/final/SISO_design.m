close all
% States: omega, z1_dot, z1, z2_dot, z2
% Inputs: beta, tau(generator torque), wind_speed
s = tf('s');
mimo = ss(A, B, C, D);
mimo.u = {'Beta', 'Tau', 'Wind Speed'};
mimo.y = {'omega', 'z'};

siso = tf(mimo);
siso = minreal(-siso(1,1));

bode(siso)
grid
figure
pzplot(siso)
grid
figure
%margin(siso)
%controlSystemDesigner('bode', siso)
%pzmap(mimo)

%% Design 1 %%
C1 = 0.050425*(s+0.0005424)/(s*(s+0.2248)*(s+0.0005843));
L1 = siso*C1;
I1 = eye(size(L1));
S1 = feedback(I1,L1); 
T1 = I1 - S1;
cloesed_loop1 = stepinfo(T1);
sensitivity1 = stepinfo(S1);
loop1 = getGainCrossover(L1, 1);
bode(S1)
grid
figure

%% Design 2 %%
C2 = 0.15299 * (s^2 + 0.01552*s + 0.03945)/(s*(s+0.1986)^2);
L2 = siso*C2;
I2 = eye(size(L2));
S2 = feedback(I2,L2); 
T2 = I2 - S2;
bode(S2)
grid
figure
cloesed_loop2 = stepinfo(T2);
sensitivity2 = stepinfo(S2);
loop2 = getGainCrossover(L2, 1)

%% Simulation and Comparison %%
step(S1, S2)
legend(["Design 1", "Design 2"])
figure
step(T1, T2)
legend(["Design 1", "Design 2"])
figure
step(S1)
figure
