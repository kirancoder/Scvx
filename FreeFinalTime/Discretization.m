classdef Discretization < handle
    % Discretization Handles First-Order Hold (FOH) discretization using ode45
    % Free-Final-Time Implementation.
    
    properties
        K
        m
        n_x
        n_u
        
        A_bar
        B_bar
        C_bar
        S_bar  % NEW: Time-dilation matrix
        z_bar
        
        % Indices for flattened vectors
        x_ind
        A_bar_ind
        B_bar_ind
        C_bar_ind
        S_bar_ind % NEW
        z_bar_ind
        
        f_func
        A_func
        B_func
        
        V0
        dt
    end
    
    methods
        function obj = Discretization(m, K)
            obj.K = K;
            obj.m = m;
            obj.n_x = m.n_x;
            obj.n_u = m.n_u;
            
            obj.A_bar = zeros(m.n_x * m.n_x, K - 1);
            obj.B_bar = zeros(m.n_x * m.n_u, K - 1);
            obj.C_bar = zeros(m.n_x * m.n_u, K - 1);
            obj.S_bar = zeros(m.n_x, K - 1);
            obj.z_bar = zeros(m.n_x, K - 1);
            
            % 1-based vector indices for flat matrices in MATLAB
            x_end = m.n_x;
            A_bar_end = m.n_x * (1 + m.n_x);
            B_bar_end = m.n_x * (1 + m.n_x + m.n_u);
            C_bar_end = m.n_x * (1 + m.n_x + m.n_u + m.n_u);
            S_bar_end = m.n_x * (1 + m.n_x + m.n_u + m.n_u + 1);
            z_bar_end = m.n_x * (1 + m.n_x + m.n_u + m.n_u + 2);
            
            obj.x_ind = 1 : x_end;
            obj.A_bar_ind = (x_end + 1) : A_bar_end;
            obj.B_bar_ind = (A_bar_end + 1) : B_bar_end;
            obj.C_bar_ind = (B_bar_end + 1) : C_bar_end;
            obj.S_bar_ind = (C_bar_end + 1) : S_bar_end;
            obj.z_bar_ind = (S_bar_end + 1) : z_bar_end;
            
            [obj.f_func, obj.A_func, obj.B_func] = m.get_equations();
            
            % Integration initial condition
            obj.V0 = zeros(z_bar_end, 1);
            obj.V0(obj.A_bar_ind) = reshape(eye(m.n_x), [], 1);
            
            obj.dt = 1.0 / (K - 1); % Normalized dt (real time is dt * sigma)
        end
        
        function [A_bar_out, B_bar_out, C_bar_out, S_bar_out, z_bar_out] = calculate_discretization(obj, X, U, sigma)
            % Calculate discretization for given states, inputs and total time
            
            options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
            
            for k = 1:(obj.K - 1)
                obj.V0(obj.x_ind) = X(:, k);
                
                u_t0 = U(:, k);
                u_t1 = U(:, k + 1);
                
                % Integrate using ode45
                [~, V_out] = ode45(@(t, V) obj.ode_dVdt(t, V, u_t0, u_t1, sigma), [0, obj.dt], obj.V0, options);
                V_final = V_out(end, :)';
                
                % Extract matrices
                Phi = reshape(V_final(obj.A_bar_ind), obj.n_x, obj.n_x);
                obj.A_bar(:, k) = reshape(Phi, [], 1);
                
                B_mat = reshape(V_final(obj.B_bar_ind), obj.n_x, obj.n_u);
                obj.B_bar(:, k) = reshape(Phi * B_mat, [], 1);
                
                C_mat = reshape(V_final(obj.C_bar_ind), obj.n_x, obj.n_u);
                obj.C_bar(:, k) = reshape(Phi * C_mat, [], 1);
                
                obj.S_bar(:, k) = Phi * V_final(obj.S_bar_ind);
                obj.z_bar(:, k) = Phi * V_final(obj.z_bar_ind);
            end
            
            A_bar_out = obj.A_bar;
            B_bar_out = obj.B_bar;
            C_bar_out = obj.C_bar;
            S_bar_out = obj.S_bar;
            z_bar_out = obj.z_bar;
        end
        
        function dVdt = ode_dVdt(obj, t, V, u_t0, u_t1, sigma)
            % ODE function to compute dVdt
            alpha = (obj.dt - t) / obj.dt;
            beta = t / obj.dt;
            
            x = V(obj.x_ind);
            u = u_t0 + (t / obj.dt) * (u_t1 - u_t0);
            
            Phi_A_xi = inv(reshape(V(obj.A_bar_ind), obj.n_x, obj.n_x));
            
            A_subs = sigma * obj.A_func(x, u);
            B_subs = sigma * obj.B_func(x, u);
            f_subs = obj.f_func(x, u);
            
            dVdt = zeros(size(V));
            dVdt(obj.x_ind) = sigma * f_subs;
            
            % Dynamics of the State Transition Matrix
            Phi_dot = A_subs * reshape(V(obj.A_bar_ind), obj.n_x, obj.n_x);
            dVdt(obj.A_bar_ind) = reshape(Phi_dot, [], 1);
            
            dVdt(obj.B_bar_ind) = reshape(Phi_A_xi * B_subs, [], 1) * alpha;
            dVdt(obj.C_bar_ind) = reshape(Phi_A_xi * B_subs, [], 1) * beta;
            
            dVdt(obj.S_bar_ind) = Phi_A_xi * f_subs;
            
            z_t = -A_subs * x - B_subs * u;
            dVdt(obj.z_bar_ind) = Phi_A_xi * z_t;
        end
        
        function X_nl = integrate_nonlinear_piecewise(obj, X_l, U, sigma)
            % Piecewise integration to verify accuracy of linearization
            X_nl = zeros(size(X_l));
            X_nl(:, 1) = X_l(:, 1);
            
            options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
            
            for k = 1:(obj.K - 1)
                u_t0 = U(:, k);
                u_t1 = U(:, k + 1);
                
                % Note: True integration time is scaled by sigma
                [~, x_out] = ode45(@(t, x) obj.dx(t, x, u_t0, u_t1, sigma), [0, obj.dt * sigma], X_l(:, k), options);
                X_nl(:, k + 1) = x_out(end, :)';
            end
        end
        
        function X_nl = integrate_nonlinear_full(obj, x0, U, sigma)
            % Simulate nonlinear behavior given an initial state and input over time
            X_nl = zeros(length(x0), obj.K);
            X_nl(:, 1) = x0;
            
            options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
            
            for k = 1:(obj.K - 1)
                u_t0 = U(:, k);
                u_t1 = U(:, k + 1);
                
                [~, x_out] = ode45(@(t, x) obj.dx(t, x, u_t0, u_t1, sigma), [0, obj.dt * sigma], X_nl(:, k), options);
                X_nl(:, k + 1) = x_out(end, :)';
            end
        end
        
        function dx_dt = dx(obj, t, x, u_t0, u_t1, sigma)
            % Simple nonlinear dynamics evaluation over real time
            u = u_t0 + (t / (obj.dt * sigma)) * (u_t1 - u_t0);
            dx_dt = obj.f_func(x, u);
        end
    end
end