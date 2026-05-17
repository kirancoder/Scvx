classdef RocketLanding3D < handle
    % RocketLanding3D A 6 degree of freedom rocket landing problem.
    
    properties
        n_x = 14;
        n_u = 3;
        K % Stored to use in constraints
        
        % Mass
        m_wet = 30000.0;
        m_dry = 22000.0;
        
        t_f_guess = 15.0;
        
        % State constraints
        r_I_init = [0.0; 200.0; 200.0];
        v_I_init = [-50.0; -100.0; -50.0];
        q_B_I_init;
        w_B_init = [0.0; 0.0; 0.0];
        
        r_I_final = [0.0; 0.0; 0.0];
        v_I_final = [0.0; 0.0; -5.0];
        q_B_I_final;
        w_B_final = [0.0; 0.0; 0.0];
        
        w_B_max;
        
        % Angles
        max_gimbal = 7;
        max_angle = 70;
        glideslope_angle = 20;
        
        tan_delta_max;
        cos_theta_max;
        tan_gamma_gs;
        
        % Thrust limits
        T_max = 800000.0;
        T_min;
        
        % Angular moment of inertia
        J_B = diag([4000000.0, 4000000.0, 100000.0]);
        
        % Gravity
        g_I = [0.0; 0.0; -9.81];
        
        % Fuel consumption
        alpha_m = 1 / (282 * 9.81);
        
        % Vector from thrust point to CoM
        r_T_B = [0.0; 0.0; -14.0];
        
        r_scale;
        m_scale;
        
        x_init;
        x_final;
    end
    
    methods
        function obj = RocketLanding3D(K_nodes)
            obj.K = K_nodes;
            obj.T_min = obj.T_max * 0.4;
            obj.w_B_max = deg2rad(90);
            
            obj.tan_delta_max = tan(deg2rad(obj.max_gimbal));
            obj.cos_theta_max = cos(deg2rad(obj.max_angle));
            obj.tan_gamma_gs = tan(deg2rad(obj.glideslope_angle));
            
            % Generate random initial state per the Python script
            % obj.r_I_init(3) = 500;
            % obj.r_I_init(1:2) = -300 + 600 * rand(2, 1);
            % obj.v_I_init(3) = -100 + 40 * rand(1);
            % obj.v_I_init(1:2) = (-0.5 + 0.3 * rand(2, 1)) .* obj.r_I_init(1:2);
            % 
            % % Initialize quaternions (euler_to_quat helper logic included here)
            % eul_init = deg2rad([-30 + 60*rand(1), -30 + 60*rand(1), 0]);
            % obj.q_B_I_init = obj.euler_to_quat(eul_init);
            % obj.q_B_I_final = obj.euler_to_quat([0, 0, 0]);
            % 
            % obj.w_B_init = deg2rad([-20 + 40*rand(1); -20 + 40*rand(1); 0]);

            % --- Fixed Initial State ---
            % Position: [East, North, Up] in meters
            obj.r_I_init = [250.0; -200.0; 500.0]; 
            
            % Velocity: [East, North, Up] in m/s
            % (Moving towards the pad and falling)
            obj.v_I_init = [-30.0; -20.0; -80.0]; 
            
            % Initial Attitude (Euler angles to Quaternions)
            % [Roll, Pitch, Yaw] -> e.g., tilted 20 degrees on X, -15 on Y
            eul_init = deg2rad([20.0; -15.0; 0.0]); 
            obj.q_B_I_init = obj.euler_to_quat(eul_init);
            obj.q_B_I_final = obj.euler_to_quat([0.0; 0.0; 0.0]);
            
            % Initial Angular Velocity: [wx, wy, wz] in rad/s
            % (Slight tumbling motion to be corrected)
            obj.w_B_init = deg2rad([-5.0; 5.0; 0.0]); 
            % ---------------------------
            
            obj.x_init = [obj.m_wet; obj.r_I_init; obj.v_I_init; obj.q_B_I_init; obj.w_B_init];
            obj.x_final = [obj.m_dry; obj.r_I_final; obj.v_I_final; obj.q_B_I_final; obj.w_B_final];
            
            obj.r_scale = norm(obj.r_I_init);
            obj.m_scale = obj.m_wet;
        end
        
        function nondimensionalize(obj)
            obj.alpha_m = obj.alpha_m * obj.r_scale;
            obj.r_T_B = obj.r_T_B / obj.r_scale;
            obj.g_I = obj.g_I / obj.r_scale;
            obj.J_B = obj.J_B / (obj.m_scale * obj.r_scale^2);
            
            obj.x_init = obj.x_nondim(obj.x_init);
            obj.x_final = obj.x_nondim(obj.x_final);
            
            obj.T_max = obj.u_nondim(obj.T_max);
            obj.T_min = obj.u_nondim(obj.T_min);
            
            obj.m_wet = obj.m_wet / obj.m_scale;
            obj.m_dry = obj.m_dry / obj.m_scale;
        end
        
        function x = x_nondim(obj, x)
            x(1) = x(1) / obj.m_scale;
            x(2:4) = x(2:4) / obj.r_scale;
            x(5:7) = x(5:7) / obj.r_scale;
        end
        
        function u = u_nondim(obj, u)
            u = u / (obj.m_scale * obj.r_scale);
        end
        
        function redimensionalize(obj)
            obj.alpha_m = obj.alpha_m / obj.r_scale;
            obj.r_T_B = obj.r_T_B * obj.r_scale;
            obj.g_I = obj.g_I * obj.r_scale;
            obj.J_B = obj.J_B * (obj.m_scale * obj.r_scale^2);
            
            obj.T_max = obj.u_redim(obj.T_max);
            obj.T_min = obj.u_redim(obj.T_min);
            
            obj.m_wet = obj.m_wet * obj.m_scale;
            obj.m_dry = obj.m_dry * obj.m_scale;
            
            obj.x_init = obj.x_redim(obj.x_init);
            obj.x_final = obj.x_redim(obj.x_final);
        end
        
        function x = x_redim(obj, x)
            x(1, :) = x(1, :) * obj.m_scale;
            x(2:4, :) = x(2:4, :) * obj.r_scale;
            x(5:7, :) = x(5:7, :) * obj.r_scale;
        end
        
        function u = u_redim(obj, u)
            u = u * (obj.m_scale * obj.r_scale);
        end
        
        function [f_func, A_func, B_func] = get_equations(obj)
            % 14 states, 3 inputs
            x_sym = sym('x_sym', [14 1], 'real');
            u_sym = sym('u_sym', [3 1], 'real');
            
            mass = x_sym(1);
            r_vec = x_sym(2:4);
            v_vec = x_sym(5:7);
            q_vec = x_sym(8:11);
            w_vec = x_sym(12:14);
            
            C_B_I = obj.dir_cosine(q_vec);
            C_I_B = C_B_I';
            
            f_sym = sym(zeros(14, 1));
            
            % Dynamics exactly mirroring Python's derivation
            f_sym(1) = -obj.alpha_m * sqrt(u_sym(1)^2 + u_sym(2)^2 + u_sym(3)^2); % norm(u)
            f_sym(2:4) = v_vec;
            f_sym(5:7) = (1 / mass) * C_I_B * u_sym + obj.g_I;
            f_sym(8:11) = 0.5 * obj.omega(w_vec) * q_vec;
            f_sym(12:14) = obj.J_B \ (obj.skew(obj.r_T_B) * u_sym) - obj.skew(w_vec) * w_vec;
            
            A_sym = jacobian(f_sym, x_sym);
            B_sym = jacobian(f_sym, u_sym);
            
            f_func = matlabFunction(f_sym, 'Vars', {x_sym, u_sym});
            A_func = matlabFunction(A_sym, 'Vars', {x_sym, u_sym});
            B_func = matlabFunction(B_sym, 'Vars', {x_sym, u_sym});
        end
        
        function [X, U] = initialize_trajectory(obj, K)
            X = zeros(obj.n_x, K);
            U = zeros(obj.n_u, K);
            
            for k = 1:K
                alpha1 = (K - k) / K;
                alpha2 = (k - 1) / K;
                
                m_k = alpha1 * obj.x_init(1) + alpha2 * obj.x_final(1);
                r_I_k = alpha1 * obj.x_init(2:4) + alpha2 * obj.x_final(2:4);
                v_I_k = alpha1 * obj.x_init(5:7) + alpha2 * obj.x_final(5:7);
                q_B_I_k = [1; 0; 0; 0];
                w_B_k = alpha1 * obj.x_init(12:14) + alpha2 * obj.x_final(12:14);
                
                X(:, k) = [m_k; r_I_k; v_I_k; q_B_I_k; w_B_k];
                U(:, k) = ((obj.T_max - obj.T_min) / 2) * [0; 0; 1];
            end
        end
        
        function model_cost = apply_constraints(obj, X_v, U_v, U_last)
            % Applies constraints and calculates the 3D-specific objective
            
            % Boundary conditions
            X_v(1, 1) == obj.x_init(1);
            X_v(2:4, 1) == obj.x_init(2:4);
            X_v(5:7, 1) == obj.x_init(5:7);
            X_v(8:11, 1) == obj.x_init(8:11);
            X_v(12:14, 1) == obj.x_init(12:14);
            
            X_v(2:14, end) == obj.x_final(2:14);
            
            % State constraints
            X_v(1, :) >= obj.m_dry; % minimum mass
            
            % Glideslope
            for k = 1:obj.K
                norm(X_v(2:3, k)) <= X_v(4, k) / obj.tan_gamma_gs;
            end
            
            % Maximum angle limit
            for k = 1:obj.K
                norm(X_v(9:10, k)) <= sqrt((1 - obj.cos_theta_max) / 2);
            end
            
            % Max angular velocity
            for k = 1:obj.K
                norm(X_v(12:14, k)) <= obj.w_B_max;
            end
            
            % Control constraints
            for k = 1:obj.K
                norm(U_v(1:2, k)) <= obj.tan_delta_max * U_v(3, k); % gimbal angle
                norm(U_v(:, k)) <= obj.T_max; % upper thrust
            end
            
            % Linearized lower thrust constraint (Introduces s_prime slack variable)
            variable s_prime(obj.K, 1) nonnegative
            for k = 1:obj.K
                u_k_last = U_last(:, k);
                n_u = norm(u_k_last);
                if n_u > 1e-6
                    obj.T_min - (u_k_last' / n_u) * U_v(:, k) <= s_prime(k);
                else
                    obj.T_min <= s_prime(k);
                end
            end
            
            % The objective to return back to SCProblem
            model_cost = 1e5 * sum(s_prime);
        end
        
        function cost = get_linear_cost(obj)
            % Simplified for the loop tracking; usually track s_prime.value here
            cost = 0; 
        end
        
        function cost = get_nonlinear_cost(obj, X, U)
            cost = 0;
            for k = 1:obj.K
                magnitude = norm(U(:, k));
                if magnitude < obj.T_min
                    cost = cost + (obj.T_min - magnitude);
                end
            end
        end
        
        % --- Math Helper Functions ---
        function q = euler_to_quat(~, a)
            cy = cos(a(2) * 0.5); sy = sin(a(2) * 0.5);
            cr = cos(a(1) * 0.5); sr = sin(a(1) * 0.5);
            cp = cos(a(3) * 0.5); sp = sin(a(3) * 0.5);
            
            q = zeros(4, 1);
            q(1) = cy * cr * cp + sy * sr * sp; % q_w
            q(2) = cy * sr * cp - sy * cr * sp; % q_x
            q(3) = sy * cr * cp + cy * sr * sp; % q_y (Fixed typo in original python where index 3 and 2 were flipped)
            q(4) = cy * cr * sp - sy * sr * cp; % q_z
        end
        
        function S = skew(~, v)
            S = [0, -v(3), v(2);
                 v(3), 0, -v(1);
                 -v(2), v(1), 0];
        end
        
        function C = dir_cosine(~, q)
            C = [1 - 2*(q(3)^2 + q(4)^2), 2*(q(2)*q(3) + q(1)*q(4)), 2*(q(2)*q(4) - q(1)*q(3));
                 2*(q(2)*q(3) - q(1)*q(4)), 1 - 2*(q(2)^2 + q(4)^2), 2*(q(3)*q(4) + q(1)*q(2));
                 2*(q(2)*q(4) + q(1)*q(3)), 2*(q(3)*q(4) - q(1)*q(2)), 1 - 2*(q(2)^2 + q(3)^2)];
        end
        
        function O = omega(~, w)
            O = [0, -w(1), -w(2), -w(3);
                 w(1), 0, w(3), -w(2);
                 w(2), -w(3), 0, w(1);
                 w(3), w(2), -w(1), 0];
        end
    end
end