%% Parameters Initialization
% Tuning parameters for q and qd
A1 = 0.25;
A2 = 0.3;
T = 2.6;
B1 = 2*pi/(T); %2.035;

% Parameter for bag bang control
maxTorque = 50;
tolerance = 0.015;`
f1 = 0;

% Parameters for Exoskeleton System
I2 = eye(2);
B = 3.17 * 10^-4;
K = 635;
C_q_dot_q = 0;
M_q = 3.12 * 10^-2;
Cd = 15 * I2;
Kd = 13 * I2;
Kz = 25 * I2;
Kv = (10^-3) * I2;
M_q_1 = 1/M_q;

% Impedance Control Parameters
lambda1 = 1;
lambda2 = 1;
chi1 = 10;
chi2 = 12;
k_g = 0.001;

% Simulation parameters
T = 40;
dt = 0.3;
time = 0:dt:T;

% Assumption data
sigma1 = 0.05; % ||M'(q)||
sigma2 = M_q; % ||M(q)||
beta = 1; % scaling factor

A_inv = 0.5 * (sigma1 + 2 * beta * sigma2);

numSteps = length(time);

z = zeros(2, numSteps);
qr_dot = zeros(2, numSteps);
torque_e = zeros(2, numSteps);

w_s = zeros(numSteps,1);

sgn = zeros(2, numSteps);
us = zeros(2, numSteps);
qr_dot_dot = zeros(2, numSteps);
q = zeros(2, numSteps);
q_dot = zeros(2, numSteps);
q_dot_dot = zeros(2, numSteps);
qd = zeros(2, numSteps);
qd_dot = zeros(2, numSteps);
qd_dot_dot = zeros(2, numSteps);

theta = zeros(2, numSteps);

torque_err =zeros(2,numSteps);
torque_e_bang_bang = zeros(2, numSteps);
us_bang_bang = zeros(2, numSteps);

z_dot= zeros(2, numSteps);
s = zeros(numSteps,1); 
qr_dot_1 = zeros(numSteps,1);
%
torque_head_e =zeros(2, numSteps);
p_qdot_q = zeros(2,numSteps);
L_q_qdot = zeros(numSteps,1);
g_q = zeros(2,numSteps);
torque_tilde_e = zeros(2,numSteps);
torque_e_dot = zeros(2,numSteps);

% bang bang control
u_bang_bang = zeros(2,numSteps);
error_Threshold = [0.1;0.1];
error = zeros(2,numSteps);

%% Simulation
for i = 1:numSteps
    % Adding f1 and f2 to simulate back the steps in reality
     % Compute state variables
    q(:, i) = [(A1 * sin(B1 * time(i) + f1)-0.4); 
               (A1 * sin(B1 * time(i) + f1 + (pi / 2))-0.4)];               
    q_dot(:, i) = [(A1 * B1 * cos(B1 * time(i) + f1)); 
                   (A1 * B1 * cos(B1 * time(i) + f1 + (pi / 2)))];
    q_dot_dot(:, i) = [(-A1 * B1 * B1 * sin(B1 * time(i) + f1)); 
                       (-A1 * B1 * B1 * sin(B1 * time(i) + f1 + (pi / 2)))];
                   
    qd(:, i) = [(A2 * sin(B1 * time(i) + f1)-0.4); 
                (A2 * sin(B1 * time(i) + f1 + (pi / 2))-0.4)];
    qd_dot(:, i) = [(A2 * B1 * cos(B1 * time(i) + f1)); 
                    (A2 * B1 * cos(B1 * time(i) + f1 + (pi / 2)))];
    qd_dot_dot(:, i) = [(-A2 * B1 * B1 * sin(B1 * time(i) + f1)); 
                        (-A2 * B1 * B1 * sin(B1 * time(i) + f1 + (pi / 2)))];
                    
    theta(:, i) = [(0.35 * sin(B1 * time(i) + f1)-0.4); 
               (0.45 * sin(B1 * time(i) + f1 + (pi / 2))-0.4)];
    
    g_q(:,i) = 2.2 * sin(q(:, i));
    
    %Visualize back the situation when there is anomaly score make the w_s
    %into 0
    if i < (numSteps)/3 || i > 2*numSteps/3
        s(i) = random('Uniform',90,110); 
    else
        s(i) = random('Uniform',120,130);
    end
    
    w_s(i,:) = lambda1 * tanh((-s(i) / chi1) + chi2) + lambda2;

    torque_e(:, i) = w_s(i,:) * (Cd * (q_dot(:, i) - qd_dot(:, i)) + Kd * (q(:, i) - qd(:, i)));
    
    torque_e_dot(:,i) = w_s(i,:) * (Cd * (q_dot_dot(:, i) - qd_dot_dot(:, i)) + Kd * (q_dot(:, i) - qd_dot(:, i)));
    
    p_qdot_q(:,i) = A_inv*q_dot(:,i);
    
    L_q_qdot(i,:) = A_inv*M_q_1;
    
    if i >= 2
        torque_tilde_e(:,i) = torque_tilde_e(:,i-1) - L_q_qdot(i,:) * torque_tilde_e(:,i-1) * dt - torque_e_dot(:,i) * dt;
        torque_head_e(:, i) = torque_e(:, i) + torque_tilde_e(:, i);
    else
        torque_head_e(:, i) = torque_e(:, i) + torque_tilde_e(:, i);
    end

    qr_dot(:, i) = qd_dot(:, i) - pinv(Cd) * Kd * (q(:, i) - qd(:, i)) + ((1/w_s(i,:))*pinv(Cd) * torque_head_e(:, i));
 
    z(:, i) = q_dot(:, i) - qr_dot(:, i); 
           
    sgn(:, i) = sign(z(:, i));
        
    us(:,i) = -Kz * z(:, i) - torque_head_e(:,i) - k_g * sgn(:, i) + (M_q + B) * qr_dot_dot(:, i) + C_q_dot_q * qr_dot(:, i) + g_q(:,i);
    
end

%% bang bang control

%% bang bang control

max_control = zeros(2, 1);
min_control = zeros(2, 1);

% Define max and min for each leg
max_control(1) = max(us(1,:)); 
min_control(1) = min(us(1,:)); 
max_control(2) = max(us(2,:)); 
min_control(2) = min(us(2,:)); 

for i = 1:numSteps
error(:,i) = z(:,i);

    % Right leg
    if abs(error(1,i)) > error_Threshold(1,1)
    if error(1,i) > 0
    u_bang_bang(1,i) = max_control(1); 
    else
    u_bang_bang(1,i) = min_control(1); 
    end
    else
    u_bang_bang(1,i) = max_control(1) * error(1,i) / error_Threshold(1,1);
    end

    % Left leg
    if abs(error(2,i)) > error_Threshold(2,1)
    if error(2,i) > 0
    u_bang_bang(2,i) = max_control(2); 
    else
    u_bang_bang(2,i) = min_control(2); 
    end
    else
    u_bang_bang(2,i) = max_control(2) * error(2,i) / error_Threshold(2,1);
    end
end

% 
% Plotting the comparison between u_bang_bang and u
figure;
subplot(2,1,1);
plot(time, u_bang_bang(1,:), 'r', 'LineWidth', 1.5); hold on;
plot(time, torque_head_e(1,:), 'b--', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Control Input for q_1');
title('Comparison of u_bang_bang and and torque_head_e (First Control Input)');
legend('u-bangbang', 'torque_head_e');
grid on;

subplot(2,1,2);
plot(time, u_bang_bang(2,:), 'r', 'LineWidth', 1.5); hold on;
plot(time, torque_head_e(2,:), 'b--', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Control Input for q_2');
title('Comparison of u_bang_bang and torque_head_e (Second Control Input)');
legend('u-bangbang', 'torque_head_e');
grid on;