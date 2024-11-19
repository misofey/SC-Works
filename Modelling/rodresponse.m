b = 1/690;
a = 1/690;
rod = tf([b], [1, a]);
sys = ss(rod);
A = sys.A;
B = sys.B;
C = sys.C;

step(rod);