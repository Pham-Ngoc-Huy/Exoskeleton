%% Parameters Initialization
% Tuning parameters for q and qd
A1 = 0.25;
A2 = 0.3;
T = 2.6;
B1 = 2*pi/(T); %2.035;

% Parameter for bag bang control
maxTorque = 50;
tolerance = 0.015;


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
s1 =zeros(2, numSteps);
s2 =zeros(2, numSteps);
us_test = zeros(2, numSteps);
flag = zeros(2, numSteps);
us_2 = zeros(2, numSteps);
z_dot= zeros(2, numSteps);
s = zeros(numSteps,1); 
qr_dot_1 = zeros(numSteps,1);

z_after = zeros(2, numSteps);
%
torque_head_e =zeros(2, numSteps);
p_qdot_q = zeros(2,numSteps);
L_q_qdot = zeros(numSteps,1);
g_q = zeros(2,numSteps);
torque_tilde_e = zeros(2,numSteps);
torque_e_dot = zeros(2,numSteps);
%% Simulation
for i = 1:numSteps
    % Adding f1 and f2 to simulate back the steps in reality
    f1 = random("Uniform",0.02,0.05); 
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
%% Plot Results
figure;
% Joint Angle vs Desired Trajectory for Left Leg
subplot(4, 2, 1);
plot(time, q(1, :)); hold on;
plot(time, qd(1, :));
title('Joint Angle vs Desired Trajectory for Left Leg');
xlabel('Time (s)');
ylabel('Angle (rad)');
legend('q', '$q_d$', 'Interpreter', 'latex');

% Joint Angle vs Desired Trajectory for Right Leg
subplot(4, 2, 2);
plot(time, q(2, :)); hold on;
plot(time, qd(2, :));
title('Joint Angle vs Desired Trajectory for Right Leg');
xlabel('Time (s)');
ylabel('Angle (rad)');
legend('q', '$q_d$', 'Interpreter', 'latex');

% Output Control and Interaction Torque for Left Leg
subplot(4, 2, 3);
plot(time, us(1,:));hold on; 
plot(time, torque_head_e(1,:));
title('Output Control and Interaction Torque for Left Leg');
xlabel('Time (s)');
ylabel('Torque (Nm)');
legend('Control Output', '$\hat{t}_e$', 'Interpreter', 'latex');

% Output Control and Interaction Torque for Right Leg
subplot(4, 2, 4);
plot(time, us(2, :));hold on;
plot(time, torque_head_e(2, :));
title('Output Control and Interaction Torque for Right Leg');
xlabel('Time (s)');
ylabel('Torque (Nm)');
legend('Control Output', '$\hat{t}_e$', 'Interpreter', 'latex');


% % Output Control and Interaction Torque for Left Leg
% subplot(4,2,5 );
% plot(time, us_bang_bang(1, :)); hold on;
% plot(time, torque_e_bang_bang(1, :));
% title('Output Control and Interaction Torque for Left Leg - Bang Bang');
% xlabel('Time (s)');
% ylabel('Deviation (rad/s)');
% legend('Control Output', 'Torque_e');
% 
% % Output Control and Interaction Torque for Right Leg
% subplot(4,2,6);
% plot(time, us_bang_bang(2, :)); hold on;
% plot(time, torque_e_bang_bang(2, :));
% title('Output Control and Interaction Torque for Right Leg - Bang Bang');
% xlabel('Time (s)');
% ylabel('Deviation (rad/s)');
% legend('Control Output', 'Torque_e');


subplot(4,2,7);

% Plot on the left y-axis
yyaxis left
plot(time, w_s, 'k-.', 'DisplayName', 'Weight');hold on;
plot(time, s/100, 'b--', 'DisplayName', 'Score/100'); 
ylabel('Weight & Score');
ylim([0 3]);  % Adjust based on your data

% Plot on the right y-axis
yyaxis right
plot(time, z(1,:), 'r--', 'DisplayName', 'z');
ylabel('z (Nm)');
ylim([-0.3 1.5]);  % Adjust based on your data

% Set x-axis label and title
xlabel('Time Window');
title('Temporal Variation of Anomaly Score');
grid on;

% Create a legend with all items from both y-axes
legend('Weight','Score/100', 'z');

% Set the main title for the entire figure
sgtitle('Variable Impedance Control Simulation');

subplot(4,2,8);

% Plot on the left y-axis
yyaxis left
plot(time, s/100, 'b--', 'DisplayName', 'Score/100'); hold on;
plot(time, w_s, 'k-.', 'DisplayName', 'Weight');
ylabel('Weight & Score');
ylim([0 3]);  % Adjust based on your data

% Plot on the right y-axis
yyaxis right
plot(time, z(2,:), 'r--', 'DisplayName', 'z');
ylabel('z (Nm)');
ylim([-0.3 1.5]);  % Adjust based on your data

% Set x-axis label and title
xlabel('Time Window');
title('Temporal Variation of Anomaly Score');
grid on;

% Create a legend with all items from both y-axes
legend('Weight','Score/100', 'z');

% Set the main title for the entire figure
sgtitle('Variable Impedance Control Simulation');
