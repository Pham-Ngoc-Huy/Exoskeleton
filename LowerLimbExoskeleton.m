classdef LowerLimbExoskeleton
    properties
        % Parameters for Exoskeleton System
        M = 0.0312;   % Inertia (kg·m^2)
        B = 0.000317; % Actuator damping (kg·m^2)
        K = 635;      % Stiffness of SEA (N·m)
        Cd = 50;      % Desired damping
        Kd = 200;     % Desired stiffness

        % Impedance Control Parameters
        lambda1 = 1.5;
        lambda2 = 1.0;
        chi1 = 0.1;
        chi2 = 1.0;
 
        % Control gains and other parameters
        L_q_qdot = 0.5;
        Kz = 10;
        k_q = 5;
        Kv = 50;
        M_q;
        g_q = 0;
        p_qdot_q = 0;
        
        % Simulation parameters
        T = 10;
        dt = 0.01;
        time;
        numSteps;
        
        % State variables
        q;
        q_dot;
        tau_e = 0;
        qd;
        qd_dot;
        z;
        sgn;
        anomaly_score;
        y;
        ydot;
        us;
        uf;
        u;
    end
    
    methods
        function obj = LowerLimbExoskeleton()
            % Initialize time and state variables
            obj.time = 0:obj.dt:obj.T;
            obj.numSteps = length(obj.time);
            obj.q = zeros(1, obj.numSteps);
            obj.q_dot = zeros(1, obj.numSteps);
            obj.qd = sin(0.5 * obj.time);
            obj.qd_dot = 0.5 * cos(0.5 * obj.time);
            obj.z = zeros(1, obj.numSteps);
            obj.sgn = zeros(1, obj.numSteps);
            obj.anomaly_score = 0.5 * (1 + sin(0.3 * obj.time));
            obj.y = zeros(1, obj.numSteps);
            obj.ydot = zeros(1, obj.numSteps);
            obj.us = zeros(1, obj.numSteps);
            obj.uf = zeros(1, obj.numSteps);
            obj.u = zeros(1, obj.numSteps);
            obj.M_q = obj.M;
        end
        
        function obj = simulate(obj)
            % Simulation loop
            for i = 2:obj.numSteps
                % Calculate weighting function w(s)
                w_s = obj.lambda1 * tanh((-obj.anomaly_score(i) / obj.chi1) + obj.chi2) + obj.lambda2;

                % Update impedance control values
                Ca = w_s * obj.Cd;
                Ka = w_s * obj.Kd;

                % Impedance control output
                tau_imp = Ca * (obj.q_dot(i-1) - obj.qd_dot(i)) + Ka * (obj.q(i-1) - obj.qd(i));
                
                % System dynamics
                q_ddot = (tau_imp - obj.tau_e) / obj.M;
                
                % Update states
                obj.q_dot(i) = obj.q_dot(i-1) + q_ddot * obj.dt;
                obj.q(i) = obj.q(i-1) + obj.q_dot(i) * obj.dt;
                
                % Calculate impedance deviation
                obj.z(i) = obj.q_dot(i) - obj.qd_dot(i) + (1 / obj.Cd) * obj.Kd * (obj.q(i) - obj.qd(i)) + obj.tau_e * (1 / obj.Cd) * (1 / w_s);
                
                % Sign function for z
                obj.sgn(i) = sign(obj.z(i));
                
                % Example term for damping (with position and velocity dependency)
                C_qdot_q = 0.05 * obj.q_dot(i) * cos(obj.q(i));

                % Derivative of the auxiliary error term
                obj.ydot(i) = -obj.L_q_qdot * obj.y(i-1) - obj.L_q_qdot * (obj.K * (obj.qd(i) - obj.q(i)) - C_qdot_q * obj.q_dot(i) - obj.g_q + obj.p_qdot_q);

                % Interaction torque estimate
                t_head_e = obj.y(i) + obj.p_qdot_q;

                % Stabilizing control action
                obj.us(i) = -obj.Kz * obj.z(i) - t_head_e - obj.k_q * obj.sgn(i) + (obj.M_q + obj.B) * 0 + C_qdot_q * obj.q_dot(i) + obj.g_q;

                % Feedforward control
                obj.uf(i) = -obj.Kv * (obj.qd_dot(i) - obj.q_dot(i));

                % Total control input
                obj.u(i) = obj.us(i) + obj.uf(i);

                % Update auxiliary tracking variable
                obj.y(i) = obj.y(i-1) + obj.ydot(i) * obj.dt;
            end
        end
        
        function plotResults(obj)
            % Plot results
            figure;
            subplot(3, 1, 1);
            plot(obj.time, obj.q, 'b', 'LineWidth', 1.5); hold on;
            plot(obj.time, obj.qd, 'r--', 'LineWidth', 1.5);
            title('Joint Angle vs. Desired Trajectory');
            xlabel('Time (s)');
            ylabel('Angle (rad)');
            legend('Actual', 'Desired');
            
            subplot(3, 1, 2);
            plot(obj.time, obj.u, 'k', 'LineWidth', 1.5);
            title('Impedance Deviation (z)');
            xlabel('Time (s)');
            ylabel('Deviation (rad/s)');
            
            subplot(3, 1, 3);
            plot(obj.time, obj.anomaly_score, 'm', 'LineWidth', 1.5);
            title('Anomaly Score');
            xlabel('Time (s)');
            ylabel('Anomaly Score');
            
            sgtitle('Variable Impedance Control Simulation');
        end
    end
end
