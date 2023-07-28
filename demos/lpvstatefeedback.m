% see https://yalmip.github.io/example/lpvstatefeedback/
clear; clc;

Anominal = [0 1 0;0 0 1;0 0 0];
B = [0;0;1];
Q = eye(3);
R = 1;

yalmip('clear');

%% Manual implementation
fprintf('\n');
disp('**** Manual implementation ****');
A1 = Anominal; A1(1,3) = -0.1;
A2 = Anominal; A2(1,3) =  0.1;

Y = sdpvar(3,3);
L = sdpvar(1,3,'full');

F = [Y >= 0];
F = [F, [-A1*Y-B*L + (-A1*Y-B*L)' Y L';Y inv(Q) zeros(3,1);L zeros(1,3) inv(R)] >= 0];
F = [F, [-A2*Y-B*L + (-A2*Y-B*L)' Y L';Y inv(Q) zeros(3,1);L zeros(1,3) inv(R)] >= 0];
optimize(F,-trace(Y))
K = value(L)*inv(value(Y))


%% Semi-Manual implementation
fprintf('\n');
disp('**** Semi-Manual implementation ****');
sdpvar t1 t2
A = A1*t1 + A2*t2;
F = [Y >=0];
F = [F, [-A*Y-B*L + (-A*Y-B*L)' Y L';Y inv(Q) zeros(3,1);L zeros(1,3) inv(R)] >= 0];
F = [F, 0 <= [t1 t2] <= 1, t1+t2 == 1, uncertain([t1 t2])];
optimize(F,-trace(Y))
K = value(L)*inv(value(Y))

%% Fully automatic implementation
fprintf('\n');
disp('**** Fully automatic implementation ****');
alpha = sdpvar(1);
A = Anominal;
Anominal(1,3) = alpha;
alpha = sdpvar(1);
A = [0 1 alpha;0 0 1;0 0 0];
F = [Y >=0];
F = [F, [-A*Y-B*L + (-A*Y-B*L)' Y L';Y inv(Q) zeros(3,1);L zeros(1,3) inv(R)] >= 0];
F = [F, -0.1 <= alpha <= 0.1, uncertain(alpha)];
optimize(F,-trace(Y))
K = value(L)*inv(value(Y))

%% Gain scheduling control
fprintf('\n');
disp('**** Gain scheduling control ****');
disp('**** Parameterized feedback matrix ****');
L0 = sdpvar(1,3);
L1 = sdpvar(1,3);
L = L0 + alpha*L1;
F = [Y >=0];
F = [F, [-A*Y-B*L + (-A*Y-B*L)' Y L';Y inv(Q) zeros(3,1);L zeros(1,3) inv(R)] >= 0];
F = [F, -0.1 <= alpha <= 0.1, uncertain(alpha)];
optimize(F,-trace(Y))

fprintf('\n');
disp('**** Parameterized feedback matrix and Lyapunov matrix ****');
Y0 = sdpvar(3,3);
Y1 = sdpvar(3,3);
Y = Y0 + alpha*Y1;

F = [Y >= 0];
F = [F, [-A*Y-B*L + (-A*Y-B*L)' Y L';Y inv(Q) zeros(3,1);L zeros(1,3) inv(R)] >= 0]
F = [F, -0.1 <= alpha <= 0.1, uncertain(alpha)];
optimize(F,-trace(Y))


