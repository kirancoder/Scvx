classdef SCProblem < handle
    
    properties
        m       
        K       
        
        % Solution Variables
        X
        U
        nu
        
        % Stored objective value for trust-region updates
        model_linear_cost
        
        % Problem Parameters
        A_bar
        B_bar
        C_bar
        z_bar
        
        X_last
        U_last
        
        weight_nu
        tr_radius
        
        verbose_solver = false;
    end
    
    methods
        function obj = SCProblem(m, K)
            obj.m = m;
            obj.K = K;
        end
        
        function set_parameters(obj, varargin)
            for i = 1:2:length(varargin)
                prop_name = varargin{i};
                prop_val = varargin{i+1};
                if isprop(obj, prop_name)
                    obj.(prop_name) = prop_val;
                else
                    fprintf('Warning: Parameter ''%s'' does not exist.\n', prop_name);
                end
            end
        end
        
        function val = get_variable(obj, name)
            if isprop(obj, name)
                val = obj.(name);
            else
                fprintf('Warning: Variable ''%s'' does not exist.\n', name);
                val = [];
            end
        end
        
        function error_flag = solve(obj)
            error_flag = false;
            cvx_clear; 
            
            K_steps = obj.K;
            nx = obj.m.n_x;
            nu_dim = obj.m.n_u;
            tr_rad = obj.tr_radius;
            w_nu = obj.weight_nu;
            
            A_b = obj.A_bar;
            B_b = obj.B_bar;
            C_b = obj.C_bar;
            z_b = obj.z_bar;
            
            x_last = obj.X_last;
            u_last = obj.U_last;
            
            % Model properties
            x_i = obj.m.x_init;
            x_f = obj.m.x_final;
            
            try
                if obj.verbose_solver
                    cvx_begin
                else
                    cvx_begin quiet
                end
                
                variable X_v(nx, K_steps)
                variable U_v(nu_dim, K_steps)
                variable nu_v(nx, K_steps - 1)
                
                dx = X_v - x_last;
                du = U_v - u_last;
                
                % Add 3D-specific slack variable if using 3D model
                is_3d = isa(obj.m, 'RocketLanding3D');
                if is_3d
                    variable s_prime(K_steps, 1) nonnegative
                end
                
                subject to
                    % 1. Universal Dynamics Constraints
                    for k = 1:(K_steps - 1)
                        A_mat = reshape(A_b(:, k), nx, nx);
                        B_mat = reshape(B_b(:, k), nx, nu_dim);
                        C_mat = reshape(C_b(:, k), nx, nu_dim);
                        
                        X_v(:, k+1) == A_mat * X_v(:, k) ...
                                     + B_mat * U_v(:, k) ...
                                     + C_mat * U_v(:, k+1) ...
                                     + z_b(:, k) ...
                                     + nu_v(:, k);
                    end
                    
                    % 2. Universal Trust Region Constraint
                    norm(dx(:), 1) + norm(du(:), 1) <= tr_rad;
                    

                    if isa(obj.m, 'RocketLanding2D')
                        X_v(1:2, 1) == x_i(1:2);
                        X_v(3:4, 1) == x_i(3:4);
                        X_v(5, 1)   == x_i(5);
                        X_v(6, 1)   == x_i(6);
                        X_v(:, end) == x_f;
                        
                        abs(X_v(5, :)) <= obj.m.t_max;
                        abs(X_v(6, :)) <= obj.m.w_max;
                        X_v(2, :) >= 0;
                        
                        abs(U_v(1, :)) <= obj.m.max_gimbal;
                        U_v(2, :) >= obj.m.T_min;
                        U_v(2, :) <= obj.m.T_max;
                        
                    % ----------------------------------------------------
                    % 3B. 3D-SPECIFIC CONSTRAINTS
                    % ----------------------------------------------------
                    elseif is_3d
                        % Boundary Conditions
                        X_v(1, 1) == x_i(1);
                        X_v(2:4, 1) == x_i(2:4);
                        X_v(5:7, 1) == x_i(5:7);
                        X_v(8:11, 1) == x_i(8:11);
                        X_v(12:14, 1) == x_i(12:14);
                        
                        X_v(2:14, end) == x_f(2:14); % Final mass is free
                        
                        % State and Control Constraints per node
                        X_v(1, :) >= obj.m.m_dry;
                        
                        for k = 1:K_steps
                            % Glideslope
                            norm(X_v(2:3, k)) <= X_v(4, k) / obj.m.tan_gamma_gs;
                            
                            % Maximum angle
                            norm(X_v(9:10, k)) <= sqrt((1 - obj.m.cos_theta_max) / 2);
                            
                            % Maximum angular velocity
                            norm(X_v(12:14, k)) <= obj.m.w_B_max;
                            
                            % Gimbal angle constraint
                            norm(U_v(1:2, k)) <= obj.m.tan_delta_max * U_v(3, k);
                            
                            % Upper thrust
                            norm(U_v(:, k)) <= obj.m.T_max;
                            
                            % Linearized lower thrust constraint
                            u_k_last = u_last(:, k);
                            n_u = norm(u_k_last);
                            if n_u > 1e-6
                                obj.m.T_min - (u_k_last' / n_u) * U_v(:, k) <= s_prime(k);
                            else
                                obj.m.T_min <= s_prime(k);
                            end
                        end
                    end
                
                % Define Objective
                if is_3d
                    sc_objective = w_nu * norm(nu_v(:), 1) + 1e5 * sum(s_prime);
                else
                    sc_objective = w_nu * norm(nu_v(:), 1);
                end
                
                minimize(sc_objective)
                
                cvx_end
                
                if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
                    obj.X = X_v;
                    obj.U = U_v;
                    obj.nu = nu_v;
                    
                    if is_3d
                        obj.model_linear_cost = 1e5 * sum(s_prime);
                    else
                        obj.model_linear_cost = 0;
                    end
                else
                    error_flag = true;
                end
                
            catch ME
                fprintf('CVX Encountered an Error: %s\n', ME.message);
                error_flag = true;
            end
        end
    end
end