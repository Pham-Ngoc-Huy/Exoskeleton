classdef bangbang
    properties
        % Parameters for Exoskeleton System
        M = 0.0312;   % Inertia (kg·m^2)
        I2= [1 0; 0 1]; % Identity matrix
        B = 3.17*10^-4; % Actuator damping (kg·m^2)
        K = 635;      % Stiffness of SEA (N·m)
        s = 0.01;
        Cd;
        Kd;
        % Impedance Control Parameters
        lambda1 = 1;
        lambda2 = 1;
        chi1 = 10;
        chi2 = 12;
        k_g = 12;

        % Control gains and other parameters
        L_q_qdot = 0.5;
        Kz;
        k_q = 5;
        Kv;
        M_q = 3.12*10^-2
        g_q;
        p_qdot_q = 0;

        % Simulation parameters
        T = 30;
        dt = 0.01;
        time;
        numSteps;

        % State variables
        q;
        theta;
        z;
        q_dot;
        torque_e;
        qd;
        qd_dot;
        qr_dot;
        Ca;
        Ka;
        y_dot;
        y;
        L_q_dot_q;
        t_head_e;
        p_q_dot_q;
        us;
        uf;
        u;
        sgn;
        qr_dot_dot_r;
        C_q_dot_q;
        
        % Bang-Bang Control Parameters
        maxTorque = 50;       % Maximum torque for bang-bang control
        tolerance = 0.05;     % Tolerance for acceptable deviation
    end

    methods
        function obj = bangbang()
            obj.time = 0:obj.dt:obj.T;
            obj.numSteps = length(obj.time);
            obj.q = (0.5*sin(0.2*obj.time+0.1) + 0.6);
            obj.q_dot = 0.5*0.2*cos(0.2*obj.time+0.1);
            obj.qd_dot = 0.5*0.3*cos(0.2*obj.time+0.1);
            obj.qd = 0.5*sin(0.2*obj.time) + 0.5;
            obj.z = zeros(1, obj.numSteps);
            obj.qr_dot = zeros(2, obj.numSteps);
            obj.torque_e = zeros(1, obj.numSteps);
            obj.Ka = zeros(1, obj.numSteps);
            obj.Ca = zeros(1, obj.numSteps);
            obj.g_q = 2.2*sin(obj.q);
            obj.sgn = zeros(1, obj.numSteps);
            obj.us = zeros(1, obj.numSteps);
            obj.uf = zeros(1, obj.numSteps);
            obj.u = zeros(1, obj.numSteps);
            obj.C_q_dot_q = zeros(1, obj.numSteps);
            obj.Cd = 15*obj.I2;
            obj.Kd = 13*obj.I2;
            obj.Kz = 25*obj.I2;
            obj.Kv = (10^-3)*obj.I2;
            obj.qr_dot_dot = -0.2*I2*(0.5*0.2 + 2*0.5*0.3)*sin(0.2*obj.time+0.1);
            obj.t_head_e = 2.5*sin(0.2*obj.time) + 2.5*cos(0.2*obj.time);
        end

        function obj = simulate(obj)
            for i = 1:obj.numSteps
                % Compute the position error
                position_error = obj.q(i) - obj.qd(i);
                
                % Bang-Bang Control Logic
                if abs(position_error) > obj.tolerance
                    % Apply max torque if outside tolerance
                    obj.u(i) = obj.maxTorque * sign(-position_error);
                else
                    % Set torque to zero within tolerance
                    obj.u(i) = 0;
                end
                
                % Update state variables if needed for visualization
                % (e.g., calculate qr_dot, z, etc. if they are being tracked)
                
                % Equation 56 (For reference)
                if obj.z(i) > 0
                    obj.sgn(i) = 1;
                elseif obj.z(i) < 0
                    obj.sgn(i) = -1;
                else
                    obj.sgn(i) = 0;
                end
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
            title('Control Output (Torque) - Bang-Bang Control');
            xlabel('Time (s)');
            ylabel('Torque (N·m)');
            legend('Control Output');

            subplot(3, 1, 3);
            plot(obj.time, obj.q - obj.qd, 'r', 'LineWidth', 1.5);
            title('Position Error');
            xlabel('Time (s)');
            ylabel('Error (rad)');
            legend('Position Error');
            
            sgtitle('Bang-Bang Control Simulation');
        end
    end
end



% đổi với bang bang - bang bang control là điều khiển bật/tắt khi thực tế
% và lí tưởng lớn hơn dung sai cho phép để ra tín hiệu cho động cơ hoạt
% động - ví dụ:
% nếu mà góc hiện tại với góc lí tưởng lệch lớn hơn dung sai cho phép là
% 0.01 thì lúc này động cơ tắt để người tự điều khiển và ngược lại thì sẽ
% bật để hỗ trợ khớp cho người dùng
% Tương tự đổi với torque (mô men xoắn) và control ouput  -> đặc biệt với
% mô men xoắn ta phải set 1 torque max để đảm bảo rằng torque sẽ không đạt
% đến mức đó do các ảnh hưởng đến động cơ và an toàn cho người dùng -> cúi
% cùng thì vẫn giống như cách hoạt động như trên 
% Cách làm này đơn giản và có thể áp dụng ở mọi mặt trận mà không lo là
% tính toán quá phức tạp - việc ta chỉ cần set 1 tolerance - ở trong bài
% này là s là được còn lại là thuật toán mình set theo logic trên sẽ làm
% việc.