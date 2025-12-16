

close all

%% Global Variables
global R1 R2 m1 m2
R1 = .25;
R2 = .381;
m1 =  .1;
m2 = 0.25;

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
    torque1 = norm(cross(r,Fg));
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
    torque1 = norm(cross(r1,Fg1));
    torque2 = norm(cross(r3,Fg2));

    sign1 = 1;
    if speed1<0
        sign1 = -1;
    end

    sign2 = 1;
    if speed2<0
        sign2 = -1;
    end
    torque3 = sign1*torque1 + sign2*(torque2);
end

%% MODEL OF AXIS
% MOTOR SHA25A-161SG-BO9B200
% L = 0.003;    % henries
% Res = 1.2;    % ohms
% Ke = 0.00428;  % V/rpm backemf  
% Kt = 0.385;   % Nm/Amps torque constant WITHOUT GEAR RATIO
% b = .001;     % viscus friction
% J = 0.037;    % kg-m^2 Moment of Inertia

L = 0.003;    % henries
Res = 1.2;    % ohms
Ke = 0.728;  % V/rpm backemf  
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

% Augmented Matricies for LQI
A_aug = [A, zeros(size(A,1),2); -C, zeros(size(C,1),2)];
B_aug = [B; zeros(2, size(B,2))];
E_aug = [E; zeros(2, size(E,2))];

%% Simulate
DT = 0.001;
T_final = 10;
t = 0;
state_0 = [0;0;0;0;0;0;0;0];
state = state_0;
states = [];
Torques = [];
Inputs = [];
T = 0:DT:T_final;
UpperCurrent = 4;
LowerCurrent = -4;
UpperVoltage = 24;
LowerVoltage = -24;
err1_integrate = 0;
err2_integrate = 0;

%% LQR FUNCTION**********************************
Q = diag([10, 100, 3, 10, 100, 3, 100000, 100000]);
R = diag([1 1]);

% Set reference points**********************************
angle1 = 50*(pi/180);
angle2 = 0*(pi/180);
theta1_ref = angle1;
theta2_ref = angle2;

% GOOD
% Q = diag([10, 100, 3, 10, 100, 3, 100000, 100000]);
% R = diag([10 10]);
Kgain = lqi(system,Q,R);

for ind = 0:DT:T_final
    
    % Calculate Error and integrate
    err1 = theta1_ref-state(1);
    err2 = theta2_ref-state(4);
    err1_integrate = err1_integrate+err1*DT;
    err2_integrate = err2_integrate+err2*DT;
    state(7) = err1_integrate;
    state(8) = err2_integrate;

    % Calculate Torque
    Torque = [CalcTorqueMotor1(state(1),state(4),state(2),state(5)); CalcTorqueMotor2(state(4),state(5))];
    Torques = [Torques, Torque];

    % Calculate Next States
    u_opt = -Kgain*state;
    u_opt_sat = [min(UpperVoltage, max(LowerVoltage, u_opt(1)));min(UpperVoltage, max(LowerVoltage, u_opt(2)))];
    Inputs = [Inputs,u_opt_sat];
    state_dot = A_aug*state+B_aug*u_opt_sat + E_aug*Torque;
    
    % Euler integration
    state = state + state_dot*DT;

    % Saturate the current
    state(3) = min(UpperCurrent, max(LowerCurrent, state(3)));
    state(6) = min(UpperCurrent, max(LowerCurrent, state(6)));
    states = [states, state];

end


%% PLOTS
figure;
subplot(2,2,1);
hold on
plot(T,states(1,:),'LineWidth',3);
plot(T,states(4,:),'LineWidth',3);
hold off
legend("Motor 1 Position", "Motor 2 Position")
ylabel("Position [Radians]")
xlim([0 T_final]);
ylim([-.5 1.75]);
ax = gca;  % Get the current axes handle
ax.FontSize = 16;

subplot(2,2,2);
hold on
plot(T,states(3,:),'LineWidth',3);
plot(T,states(6,:),'LineWidth',3);
hold off
legend("Motor 1 current", "Motor 2 current")
ylabel("Current [Amps]")
ax = gca;  % Get the current axes handle
ax.FontSize = 16;

subplot(2,2,3)
hold on
plot(T,Inputs(1,:),'LineWidth',3);
plot(T,Inputs(2,:),'LineWidth',3);
hold off
legend("Motor 1 Input", "Motor 2 Input")
ylabel("Voltage [Volts]")
ax = gca;  % Get the current axes handle
ax.FontSize = 16;