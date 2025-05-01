%% Clear all
clear all; close all; clc;

%% Init system
J1 = 0.0025;         %kgm^2
J2 = 0.0018;         %kgm^2
J3 = 0.0018;         %kgm^2
k1 = 2.7;            %Nm/rad
k2 = 2.6;            %Nm/rad
b1 = 0.0029;         %Nms/rad
b2 = 0.0002;         %Nms/rad
b3 = 0.00015;        %Nms/rad

A = [   0      1       0         0     0     0;
     -k1/J1 -b1/J1   k1/J1       0     0     0;
        0      0       0         1     0     0;
      k1/J2    0  -(k1+k2)/J2 -b2/J2 k2/J2   0;
        0      0       0         0     0     1;
        0      0     k2/J3       0  -k2/J3 -b3/J3];

B = [ 0  ;
     1/J1;
      0  ;
      0  ;
      0  ;
      0  ];

Ex = [  0  ;
      -1/J1;
        0  ;
        0  ;
        0  ;
        0  ];

C = [0 0 0 0 1 0];

D = zeros(1,1);

Ey = zeros(1,1);

%% Question 9

s = tf('s');
G = minreal(C*inv(s*eye(size(A,1)) - A)*B + D);
Gd = minreal(C*inv(s*eye(size(A,1)) - A)*Ex + Ey);

%Weight matrices
% % % % % % % % % %shape S to be small inside the control bandwidth
% % % % % % % % % M = 2; %M must be larger than 1
% % % % % % % % % wb = 10; % bandwidth of the CL. wb~settling time (wb = 4/settling time)
% % % % % % % % % A_ = 0.001; %low steady state error
% % % % % % % % % WS = (s/M+wb)/(s+wb*A_);
% % % % % % % % % WS = WS^4; %high order needed to have large dc gain
% % % % % % % % % 
% % % % % % % % % %frequencies where the bode plot breaks
% % % % % % % % % % z1 = 0;
% % % % % % % % % % z2 = 0;
% % % % % % % % % % z3 = 0.1;
% % % % % % % % % % p1 = 1000;
% % % % % % % % % % p2 = 1000;
% % % % % % % % % % p3 = 1000;
% % % % % % % % % % WKS = 1*(s+z1)*(s+z2)*(s+z3)/(s+p1)/(s+p2)/(s+p3); 
% % % % % % % % % % %(s+w2/M2)/(es+w2)
% % % % % % % % % M = 2; %M must be larger than 1
% % % % % % % % % wb = 5; % bandwidth of the CL. wb~settling time (wb = 4/settling time)
% % % % % % % % % A_ = 0.001; %low steady state error
% % % % % % % % % WKS = (s+wb/M)/(A_*s+wb);


% % M = 2;
% % wb = 6;
% % A_ = 0.001;
% % WS = (s/M+wb)/(s+wb*A_);
% % WS = WS^4;
% % M = 10;
% % wb = 500;
% % A_ = 0.1;
% % WKS = (s+wb/M)/(A_*s+wb);

% M = 1.7;
% wb = 6;
% A_ = 0.001;
% WS = (s/M+wb)/(s+wb*A_);
% WS = WS^3;
% M = 50;
% wb = 10000;
% A_ = 0.01;
% WKS = (s+wb/M)/(A_*s+wb);

% M = 1.6;
% wb = 6;
% A_ = 0.001;
% WS = (s/M+wb)/(s+wb*A_);
% WS = WS^3;
% M = 1000;
% wb = 5000;
% A_ = 0.1;
% WKS = (s+wb/M)/(A_*s+wb);
% WKS = WKS^2;

M = 2;
wb = 2.1;
A_ = 0.01;
WS = (s/M+wb)/(s+wb*A_);
WS = WS^3;

M = 7;
wb = 7;
A_ = 0.07;
WKS = (s+wb/M)/(A_*s+wb);
WKS = WKS^2;

log_min = -2;
log_max = 3;

%Formulate P matrix

A = A;
B2 = B;
B1 = zeros(6,2);
C2 = -C;
C1 = [-WS*C; zeros(1,6)];
D11 = [WS -WS*Gd; 0 0];
D12 = [0; WKS];
D21 = [1 -Gd];
D22 = 0;

PA = [A B1 B2;
       C1 D11 D12;
       C2 D21 D22];

%Assumptions
if (rank(ctrb(A,B2))==size(A,1))
    disp("System is controllable -> stabilizable as well.")
end
if (rank(obsv(A,C2))==size(A,1))
    disp("System is observable -> detectable as well.")
end

%just checked them
D12;
D21;
disp("D12 has full rank")
disp("D21 has full rank")

[A-s*eye(size(A)) B2; C1 D12];
[A-s*eye(size(A)) B1; C2 D21];

D11;
D22;

P = [ WS -WS*Gd -WS*G ;
      0     0    WKS  ;
      1    -Gd   -G   ];

[K_inf, CL, Gam, ~] = mixsyn(G, WS, WKS, []);
[K_inf, CL, Gam, ~] = mixsyn(G, WS, WKS, [], [Gam*1.02 Gam*1.02]); %adjust gamma to move away from optimal solution to avoid high gains

S = 1/(1+K_inf*G);
T = K_inf*G*S;

figure(1)
bode(WS, logspace(log_min, log_max, 50)); hold on;
bode(WKS, logspace(log_min, log_max, 50)); hold off;
grid on;
% step(T);
% subplot(1,4,1);
% bode(T, logspace(log_min, log_max, 50)); hold on;
% bode(S, logspace(log_min, log_max, 50)); hold off;
% legend({"T", "S"});
% grid on;
% subplot(1,4,2);
[mag, phase, omega] = bode(K_inf, logspace(log_min, log_max, 50));
mag = squeeze(mag);  % Convert magnitude matrix to a vector
omega = squeeze(omega);  % Convert omega matrix to a vector
figure(2);
semilogx(omega, 20*log10(mag)); hold on;
plot([1 1], [min(20*log10(mag)) max(20*log10(mag))], 'r', 'LineWidth', 2); hold on;
plot([min(omega) max(omega)], [0 0], 'b', 'LineWidth', 2); hold off;
legend({"K_inf"});
grid on;
% subplot(1,4,3);
figure(3);
bode(Gam/WS, logspace(log_min, log_max, 50)); hold on;
bode(S, logspace(log_min, log_max, 50)); hold off;
legend({"1/Ws", "S"});
grid on;
figure(4);
step(K_inf*S-K_inf*Gd*S);
grid on;

% subplot(1,4,4);
% bode(1/WKS, logspace(log_min, log_max, 50)); hold on;
% bode(K_inf*S, logspace(log_min, log_max, 50)); hold off;
% legend({"1/Wks", "t"});
% grid on;
Gam