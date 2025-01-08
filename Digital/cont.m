s = tf('s');
G = 1 / (s * (s + 2) * (s*s + 100*s + 2600));


%% System definition
Kp = 2.94e4;
C_manual = Kp*10*(s+2)/(s+20);

plant = G;
plant.u = 'u';
plant.y = 'y';

controller = C_manual;
controller.u = 'e';
controller.y = 'u';

err = sumblk('e = r - y');

combined = connect(controller, plant, err, 'r', 'y');

% pole(combined)
% zero(combined)
% step(combined)
damp(combined)
stepinfo(combined, 'SettlingTimeThreshold', 0.01)

%% Disturbance rejection
K2 = 1100000;
integrated = K2*(s+2)*(s+2)/(s+20)/s;

controller_disturbed = integrated;
controller_disturbed.y = 'uC';
controller_disturbed.u = 'e';
dist = sumblk('u = du + uC');
err = sumblk('e = r - y');
disturbed = connect(controller_disturbed, plant, err, dist, 'du', 'y');
step(disturbed)
stepinfo(disturbed, 'SettlingTimeThreshold', 0.01)

disturbed_complete = connect(controller_disturbed, plant, err, dist, {'r', 'du'}, 'y');
step(disturbed_regular(1))


%% controller iterations

C_new = 17717 * (s+5) / (s+16);
C_new.u = 'e';
C_new.y = 'u';
new_system = connect(C_new, plant, err, 'u', {'y', 'u'});
stepinfo(new_system(1))
step(new_system(2))