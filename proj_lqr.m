

close all

%% Global Variables
global R1 R2 m1 m2
R1 = .5;
R2 = .75;
m1 = .5;
m2 = 1.2;

%% Torque Calculation Functions

% Calculate Torque 1
function torque1=CalcTorqueMotor2(theta2,speed)
    global R2;
    global m2;
    r = [R2/2*cos(theta2), R2/2*sin(theta2), 0];
    Fg = [0, -m2*9.8, 0];
    sign = 1;
    if speed<0
        sign = -1;
    end
    torque1 = sign*norm(cross(r,Fg));
end

% Calculate Torque 1
function torque3=CalcTorqueMotor1(theta1,theta2,speed1,speed2)
    global R1;
    global R2;
    global m1;
    global m2;
    r1 = [R1/2*cos(theta1), R1/2*sin(theta1), 0];
    r2 = [R2/2*cos(theta2), R2/2*sin(theta2), 0];
    r3 = 2*r1 + r2;
    Fg1 = [0, -m1*9.8, 0];
    Fg2 = [0, -m2*9.8, 0];
    torque1 = cross(r1,Fg1);
    torque2 = cross(r3,Fg2);
    sign = 1;
    if speed1<0
        sign = -1;
    end
    sign2 = 1;
    if speed2<0
        sign2 = -1;
    end
    torque3 = norm(sign*torque1 + sign*torque2);
end

%% MODEL OF AXIS
% MOTOR SHA25A-161SG-BO9B200
% L = 0.003;    % henries
% Res = 1.2;    % ohms
% Ke = 0.00428;  % V/rad/s backemf  
% Kt = 0.385;   % Nm/Amps torque constant WITHOUT GEAR RATIO
% b = .001;     % viscus friction
% J = 0.037;    % kg-m^2 Moment of Inertia

L = 0.003;    % henries
Res = 1.2;    % ohms
Ke = 0.728;  % V/rad/s backemf  
Kt = 10.385;   % Nm/Amps torque constant WITHOUT GEAR RATIO
b = .001;     % viscus friction
J = 0.037;    % kg-m^2 Moment of Inertia
% Gear ratio 1:161

%% Setup State Space model
A = [0 1 0 0 0 0;
     0 -b/J Kt/J 0 0 0;
     0 -Ke/L -Res/L 0 0 0;
     0 0 0 0 1 0;
     0 0 0 0 -b/J Kt/J;
     0 0 0 0 -Ke/L -Res/L];

B = [0 0;0 0;1/L 0;0 0;0 0;0 1/L];
E = [0 0;-1/J 0;0 0;0 0;0 -1/J;0 0];
C = [1 0 0 0 0 0;
     0 0 0 1 0 0];
D = [0 0;0 0];
system = ss(A,B,C,D);

%% LQR FUNCTION
Q = diag([1000, 10, 10, 1000, 10, 10]);

R = [3 0;
     0 3];
Kgain = lqr(system,Q,R);

%% Simulate
angle = 70*(pi/180);
theta1_ref = angle;
theta2_ref = 0;
DT = 0.001;
T_final = 50;
t = 0;
state_0 = [0;0;0;0;0;0];
state_old = state_0;
states = [];
Torques = [];
Inputs = [];
T = 0:DT:T_final;
UpperCurrent = 4;
LowerCurrent = -4;
UpperVoltage = 24;
LowerVoltage = -24;
Uss1 = (Res/Kt)*CalcTorqueMotor1(theta1_ref,theta2_ref,0,0);
Uss2 = (Res/Kt)*CalcTorqueMotor2(theta2_ref,0);
state_ref = [theta1_ref;0;Uss1/Res;theta2_ref;0;Uss2/Res];

for ind = 0:DT:T_final

    % Solve for torque
    Torque = [CalcTorqueMotor1(state(1),state(4),state(2),state(5)); CalcTorqueMotor2(state(4),state(5))];
    Torques = [Torques, Torque];
    
    % Solve for state_dot
    u_opt = [Uss1;Uss2]-Kgain*(state_old-state_ref);


    u_opt_saturate(1) = min(UpperVoltage, max(LowerVoltage, u_opt(1)));
    u_opt_saturate(2) = min(UpperVoltage, max(LowerVoltage, u_opt(2)));
    %state_dot = A*(state_old-state_ref)+B*([Uss1;Uss2]-Kgain*(state_old-state_ref)) + E*Torque;
    state_dot = A*(state_old-state_ref)+B*(transpose(u_opt_saturate)) + E*Torque;
    Inputs = [Inputs, transpose(u_opt_saturate)];
    
    % Euler Integration
    state = state_old + state_dot*DT;

    % Saturate Current
    state(3) = min(UpperCurrent, max(LowerCurrent, state(3)));
    state(6) = min(UpperCurrent, max(LowerCurrent, state(6)));

    % Log the states
    states = [states, state];

    % Set the old states
    state_old = state;


end

%% PLOTS
figure;
hold on
plot(T,states(1,:),'LineWidth',3);
plot(T,states(4,:),'LineWidth',3);
hold off
legend("Motor 1 Position", "Motor 2 Position")
xlim([0 T_final]);
ylim([-1.75 1.75]);
ax = gca;  % Get the current axes handle
ax.FontSize = 16;

figure;
hold on
plot(T,states(3,:),'LineWidth',3);
plot(T,states(6,:),'LineWidth',3);
hold off
legend("Motor 1 current", "Motor 2 current")
ax = gca;  % Get the current axes handle
ax.FontSize = 16;

figure;
hold on
plot(T,Inputs(1,:),'LineWidth',3);
plot(T,Inputs(2,:),'LineWidth',3);
hold off
legend("Motor 1 Voltage", "Motor 2 Voltage")
ax = gca;  % Get the current axes handle
ax.FontSize = 16;