classdef RocketLanding2D < handle
    % RocketLanding2D A 2D path rocket landing problem.
    
    properties
        n_x = 6;
        n_u = 2;
        
        m = 2;
        I = 1e-2;
        
        r_init = [4.0; 4.0];
        v_init = [-2.0; -1.0];
        t_init = 0.0;
        w_init = 0.0;
        
        r_final = [0.0; 0.0];
        v_final = [0.0; 0.0];
        t_final = 0.0;
        w_final = 0.0;
        
        t_f_guess = 10; % seconds
        
        t_max;
        w_max;
        max_gimbal;
        T_max = 5;
        T_min;
        r_T = 1e-2;
        g = 1;
        
        r_scale;
        m_scale;
        
        x_init;
        x_final;
    end
    
    methods
        function obj = RocketLanding2D()
            obj.t_max = deg2rad(60);
            obj.w_max = deg2rad(60);
            obj.max_gimbal = deg2rad(7);
            obj.T_min = obj.T_max * 0.4;
            
            obj.r_scale = norm(obj.r_init);
            obj.m_scale = obj.m;
            
            obj.x_init = [obj.r_init; obj.v_init; obj.t_init; obj.w_init];
            obj.x_final = [obj.r_final; obj.v_final; obj.t_final; obj.w_final];
        end
        
        function nondimensionalize(obj)
            obj.r_init = obj.r_init / obj.r_scale;
            obj.v_init = obj.v_init / obj.r_scale;
            obj.r_T = obj.r_T / obj.r_scale;
            obj.g = obj.g / obj.r_scale;
            obj.I = obj.I / (obj.m_scale * obj.r_scale^2);
            obj.m = obj.m / obj.m_scale;
            obj.T_min = obj.T_min / (obj.m_scale * obj.r_scale);
            obj.T_max = obj.T_max / (obj.m_scale * obj.r_scale);
            
            obj.x_init = obj.x_nondim(obj.x_init);
            obj.x_final = obj.x_nondim(obj.x_final);
        end
        
        function x = x_nondim(obj, x)
            x(1:4) = x(1:4) / obj.r_scale;
        end
        
        function u = u_nondim(obj, u)
            u(2,:) = u(2,:) / (obj.m_scale * obj.r_scale);
        end
        
        function redimensionalize(obj)
            obj.r_init = obj.r_init * obj.r_scale;
            obj.v_init = obj.v_init * obj.r_scale;
            obj.r_T = obj.r_T * obj.r_scale;
            obj.g = obj.g * obj.r_scale;
            obj.I = obj.I * (obj.m_scale * obj.r_scale^2);
            obj.m = obj.m * obj.m_scale;
            obj.T_min = obj.T_min * (obj.m_scale * obj.r_scale);
            obj.T_max = obj.T_max * (obj.m_scale * obj.r_scale);
            
            obj.x_init = obj.x_redim(obj.x_init);
            obj.x_final = obj.x_redim(obj.x_final);
        end
        
        function x = x_redim(obj, x)
            x(1:4, :) = x(1:4, :) * obj.r_scale;
        end
        
        function u = u_redim(obj, u)
            u(2, :) = u(2, :) * (obj.m_scale * obj.r_scale);
        end
        
        function [f_func, A_func, B_func] = get_equations(obj)
            % FIXED: Using sym('x') instead of syms to avoid warning conflicts
            rx = sym('rx', 'real'); ry = sym('ry', 'real');
            vx = sym('vx', 'real'); vy = sym('vy', 'real');
            t_sym = sym('t_sym', 'real'); w_sym = sym('w_sym', 'real');
            gimbal = sym('gimbal', 'real'); T_sym = sym('T_sym', 'real');
            
            x_sym = [rx; ry; vx; vy; t_sym; w_sym];
            u_sym = [gimbal; T_sym];
            
            f_sym = sym(zeros(6, 1));
            f_sym(1) = x_sym(3);
            f_sym(2) = x_sym(4);
            f_sym(3) = (1 / obj.m) * sin(x_sym(5) + u_sym(1)) * u_sym(2);
            f_sym(4) = (1 / obj.m) * (cos(x_sym(5) + u_sym(1)) * u_sym(2) - obj.m * obj.g);
            f_sym(5) = x_sym(6);
            f_sym(6) = (1 / obj.I) * (-sin(u_sym(1)) * u_sym(2) * obj.r_T);
            
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
                k_py = k - 1; 
                alpha1 = (K - k_py) / K;
                alpha2 = k_py / K;
                
                X(:, k) = alpha1 * obj.x_init + alpha2 * obj.x_final;
            end
            
            U(1, :) = 0;
            U(2, :) = (obj.T_max - obj.T_min) / 2;
        end
        
        function apply_constraints(obj, X_v, U_v)
            % FIXED: Extract properties to local arrays FIRST to prevent 
            % CVX overloaded-operator dispatch errors (the 'newcnstr' error).
            xi = obj.x_init;
            xf = obj.x_final;
            tm = obj.t_max;
            wm = obj.w_max;
            mg = obj.max_gimbal;
            t_min = obj.T_min;
            t_max = obj.T_max;
            
            X_v(1:2, 1) == xi(1:2);
            X_v(3:4, 1) == xi(3:4);
            X_v(5, 1)   == xi(5);
            X_v(6, 1)   == xi(6);
            
            X_v(:, end) == xf;
            
            abs(X_v(5, :)) <= tm;
            abs(X_v(6, :)) <= wm;
            X_v(2, :) >= 0;
            
            abs(U_v(1, :)) <= mg;
            U_v(2, :) >= t_min;
            U_v(2, :) <= t_max;
        end
        
        function cost = get_linear_cost(obj)
            cost = 0;
        end
        
        function cost = get_nonlinear_cost(obj, X, U)
            cost = 0;
        end
    end
end