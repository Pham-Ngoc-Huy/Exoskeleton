classdef LowerLimbExoskeleton
    properties
        % Tunning parameter for q and qd
        A1 = 0.02;
        A2 = 0.03;
        B1 = 1;
        f1;
        f2;

        % Parameters for Exoskeleton System
        I2= eye(2); 
        B = 3.17*10^-4; 
        K = 635;      
        s = 0.001;
        C_q_dot_q = 0;
        M_q = 3.12*10^-2
        Cd;
        Kd;

        % Impedance Control Parameters
        lambda1 = 1;
        lambda2 = 1;
        chi1 = 10;
        chi2 = 12;
        k_g = 0.001;

        % Control gains and other parameters
        Kz;
        Kv;
        g_q;

        % Simulation parameters
        T = 30;
        dt = 0.1;
        time;
        numSteps;

        % State variables
        q;
        z;
        q_dot;
        torque_e;
        qd;
        qd_dot;
        qr_dot;
        Ca;
        Ka;
        t_head_e;
        p_q_dot_q;
        us;
        sgn;
        qr_dot_dot;
        w_s;
        q_dot_dot;
        qd_dot_dot;

    end

    methods
        function obj = LowerLimbExoskeleton()
            obj.f1 = rand(1);
            obj.f2 = rand(1);
            obj.time = 0:obj.dt:obj.T;
            obj.numSteps = length(obj.time);
            obj.q = [(obj.A1 * sin(obj.B1 * obj.time + obj.f1)); 
                     (obj.A1 * sin(obj.B1 * obj.time + obj.f1 + (pi / 2)))];

            obj.q_dot = [(obj.A1*obj.B1*cos(obj.B1*obj.time+ obj.f1));
                        (obj.A1*obj.B1*cos(obj.B1*obj.time+ obj.f1+ (pi/2)))];
            
            obj.q_dot_dot = [(-obj.A1*obj.B1*obj.B1*sin(obj.B1*obj.time+ obj.f1));
                             (-obj.A1*obj.B1*obj.B1*sin(obj.B1*obj.time+ obj.f1 +(pi/2)))];

            obj.qd = [(obj.A2*sin(obj.B1*obj.time + obj.f2));
                    (obj.A2*sin(obj.B1*obj.time + obj.f2 + (pi/2)))];

            obj.qd_dot = [(obj.A2*obj.B1*cos(obj.B1*obj.time + obj.f2));
                          (obj.A2*obj.B1*cos(obj.B1*obj.time + obj.f2 + (pi/2)))];
            
            obj.qd_dot_dot = [(-obj.A2*obj.B1*obj.B1*sin(obj.B1*obj.time + obj.f2));
                                -obj.A2*obj.B1*obj.B1*sin(obj.B1*obj.time + obj.f2 + (pi/2))];

            obj.g_q = 2.2 * sin(obj.q);                   
            obj.Cd = 15 * obj.I2;                          
            obj.Kd = 13 * obj.I2;                         
            obj.Kz = 25 * obj.I2;                          
            obj.Kv = (10^-3) * obj.I2;

            obj.z = zeros(2, obj.numSteps);

            obj.qr_dot = zeros(2, obj.numSteps);           
            obj.torque_e = zeros(2, obj.numSteps);         
            obj.Ka = zeros(2, 2, obj.numSteps);           
            obj.Ca = zeros(2, 2, obj.numSteps);            
            obj.sgn = zeros(2, obj.numSteps);              
            obj.us = zeros(2, obj.numSteps);               

            obj.t_head_e = obj.torque_e;   
        end

        function obj = simulate(obj)
            for i = 1:obj.numSteps
                
                obj.f1 = rand(1);
                obj.f2 = rand(1);
                
                obj.q = [(obj.A1 * sin(obj.B1 * obj.time + obj.f1)); 
                     (obj.A1 * sin(obj.B1 * obj.time + obj.f1 + (pi / 2)))];

                obj.q_dot = [(obj.A1*obj.B1*cos(obj.B1*obj.time+ obj.f1));
                            (obj.A1*obj.B1*cos(obj.B1*obj.time+ obj.f1+ (pi/2)))];

                obj.q_dot_dot = [(-obj.A1*obj.B1*obj.B1*sin(obj.B1*obj.time+ obj.f1));
                                 (-obj.A1*obj.B1*obj.B1*sin(obj.B1*obj.time+ obj.f1 +(pi/2)))];

                obj.qd = [(obj.A2*sin(obj.B1*obj.time + obj.f2));
                        (obj.A2*sin(obj.B1*obj.time + obj.f2 + (pi/2)))];

                obj.qd_dot = [(obj.A2*obj.B1*cos(obj.B1*obj.time + obj.f2));
                              (obj.A2*obj.B1*cos(obj.B1*obj.time + obj.f2 + (pi/2)))];

                obj.qd_dot_dot = [(-obj.A2*obj.B1*obj.B1*sin(obj.B1*obj.time + obj.f2));
                                    -obj.A2*obj.B1*obj.B1*sin(obj.B1*obj.time + obj.f2 + (pi/2))];
                                
                obj.g_q = 2.2 * sin(obj.q);                   

                obj.w_s = obj.lambda1 * tanh((-obj.s / obj.chi1) + obj.chi2) + obj.lambda2;
                
                obj.qr_dot_dot = obj.qd_dot_dot - inv(obj.Cd)*obj.Kd*(obj.q_dot-obj.qd_dot) + inv(obj.Cd)*obj.Cd*(obj.q_dot_dot - obj.qd_dot_dot) + obj.Kd*(obj.q_dot-obj.qd_dot);

                % Equation 42
                obj.Ca(:,:,i) = obj.w_s * obj.Cd;
                obj.Ka(:,:,i) = obj.w_s * obj.Kd;
                
                obj.torque_e(:,i) = obj.w_s*(obj.Cd * (obj.q_dot(:,i) - obj.qd_dot(:,i)) + obj.Kd * (obj.q(:,i) - obj.qd(:,i)));

                % Equation 44
                obj.qr_dot(:,i) = obj.qd_dot(:,i) - inv(obj.Cd) * obj.Kd * (obj.q(:,i) - obj.qd(:,i)) + ((1/obj.w_s) * inv(obj.Cd) * obj.torque_e(:,i));
               
                % Equation 43
                obj.z(:,i) = obj.q_dot(:,i) - obj.qr_dot(:,i);

                % Equation 56
                obj.sgn(:, i) = sign(obj.z(:, i));

                
                % Equation 55
                obj.us(:,i) = -obj.Kz*obj.z(:,i) - obj.torque_e(:,i) - obj.k_g*obj.sgn(:,i) + (obj.M_q + obj.B)*obj.qr_dot_dot(:,i) + obj.C_q_dot_q * obj.qr_dot(:,i) + obj.g_q(:,i);
               
            end
       
        end

        function plotResults(obj)
            % Plot_q_qdesired_left_leg
            figure;
            subplot(2, 2, 1);
            plot(obj.time, obj.q(1,:), 'b', 'LineWidth', 1.5); hold on;
            plot(obj.time, obj.qd(1,:), 'r--', 'LineWidth', 1.5);
            title('Joint Angle vs. Desired Trajectory for Left leg');
            xlabel('Time (s)');
            ylabel('Angle (rad)');
            legend('Actual', 'Desired');
            

            % Plot_q_qdesired_right_leg
            subplot(2, 2, 2);
            plot(obj.time, obj.q(2,:), 'b', 'LineWidth', 1.5); hold on;
            plot(obj.time, obj.qd(2,:), 'r--', 'LineWidth', 1.5);
            title('Joint Angle vs. Desired Trajectory for Right Leg');
            xlabel('Time (s)');
            ylabel('Angle (rad)');
            legend('Actual', 'Desired');

            % Plot u and t_head_e for left Leg
            subplot(2, 2, 3);
            plot(obj.time, obj.us(1,:), 'k', 'LineWidth', 1.5); hold on;
            plot(obj.time, obj.torque_e(1,:), 'b', 'LineWidth', 1.5); 
            title('Output control and Interaction Torque of Left Leg');
            xlabel('Time (s)');
            ylabel('Deviation (rad/s)')
            legend('Control Output','Torque_e');

            % Plot u and t_head_e for right Leg
            subplot(2, 2, 4);
            plot(obj.time, obj.us(2,:), 'k', 'LineWidth', 1.5); hold on;
            plot(obj.time, obj.torque_e(2,:), 'b', 'LineWidth', 1.5); 
            title('Output control and Interaction Torque of Right Leg');
            xlabel('Time (s)');
            ylabel('Deviation (rad/s)');
            legend('Control Output','Torque_e');

            sgtitle('Variable Impedance Control Simulation');
        end
    end
end

