classdef SCProblem < handle
    % SCProblem Defines a Free-Final-Time Successive Convexification problem.
    % Dynamically handles constraints for 2D and 3D Rocket Landing models.
    
    properties
        m       % Model object
        K       % Number of discretization points
        
        % Solution Variables
        X
        U
        nu
        sigma   % NEW: Free final time variable
        
        % Stored objective value for trust-region updates
        model_linear_cost
        
        % Problem Parameters (Transition Matrices)
        A_bar
        B_bar
        C_bar
        S_bar   % NEW: Time-dilation transition matrix
        z_bar
        
        % Previous Iteration States
        X_last
        U_last
        sigma_last % NEW: Previous iteration's flight time
        
        % Weights & Trust Region
        weight_sigma % NEW: Penalty weight for flight time
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
            % Assign named parameters easily
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
            % Fetch solved CVX variables
            if isprop(obj, name)
                val = obj.(name);
            else
                fprintf('Warning: Variable ''%s'' does not exist.\n', name);
                val = [];
            end
        end
        
        function error_flag = solve(obj)
            error_flag = false;
            cvx_clear; % Prevent hanging constraints from previous runs
            
            % Extract basic dimensions and scalar parameters
            K_steps = obj.K;
            nx = obj.m.n_x;
            nu_dim = obj.m.n_u;
            tr_rad = obj.tr_radius;
            w_nu = obj.weight_nu;
            w_sigma = obj.weight_sigma; % NEW
            
            % Extract matrices
            A_b = obj.A_bar;
            B_b = obj.B_bar;
            C_b = obj.C_bar;
            S_b = obj.S_bar; % NEW
            z_b = obj.z_bar;
            
            x_last = obj.X_last;
            u_last = obj.U_last;
            sig_last = obj.sigma_last; % NEW
            
            % Model boundary properties (for local parser visibility)
            x_i = obj.m.x_init;
            x_f = obj.m.x_final;
            
            try
                if obj.verbose_solver
                    cvx_begin
                else
                    cvx_begin quiet
                end
                
                % Define Decision Variables
                variable X_v(nx, K_steps)
                variable U_v(nu_dim, K_steps)
                variable nu_v(nx, K_steps - 1)
                variable sigma_v nonnegative  % NEW: Flight time variable
                
                % Trust region deltas
                dx = X_v - x_last;
                du = U_v - u_last;
                ds = sigma_v - sig_last;      % NEW: Time delta
                
                % Add 3D-specific slack variable if using 3D model
                is_3d = isa(obj.m, 'RocketLanding3D');
                if is_3d
                    variable s_prime(K_steps, 1) nonnegative
                end
                
                subject to
                    % 1. Universal Dynamics Constraints (Now includes S_bar * sigma)
                    for k = 1:(K_steps - 1)
                        A_mat = reshape(A_b(:, k), nx, nx);
                        B_mat = reshape(B_b(:, k), nx, nu_dim);
                        C_mat = reshape(C_b(:, k), nx, nu_dim);
                        S_vec = reshape(S_b(:, k), nx, 1);
                        
                        X_v(:, k+1) == A_mat * X_v(:, k) ...
                                     + B_mat * U_v(:, k) ...
                                     + C_mat * U_v(:, k+1) ...
                                     + S_vec * sigma_v ...      % NEW
                                     + z_b(:, k) ...
                                     + nu_v(:, k);
                    end
                    
                    % 2. Universal Trust Region Constraint (Includes ds)
                    norm(dx(:), 1) + norm(du(:), 1) + norm(ds, 1) <= tr_rad;
                    
                    % ----------------------------------------------------
                    % 3A. 2D-SPECIFIC CONSTRAINTS
                    % ----------------------------------------------------
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
                        X_v(1, 1) == x_i(1);
                        X_v(2:4, 1) == x_i(2:4);
                        X_v(5:7, 1) == x_i(5:7);
                        X_v(8:11, 1) == x_i(8:11);
                        X_v(12:14, 1) == x_i(12:14);
                        
                        X_v(2:14, end) == x_f(2:14); % Final mass is free
                        
                        X_v(1, :) >= obj.m.m_dry; % Minimum mass
                        
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
                
                % Define Objective (Now minimizes sigma multiplied by weight_sigma)
                if is_3d
                    sc_objective = w_sigma * sigma_v + w_nu * norm(nu_v(:), 1) + 1e5 * sum(s_prime);
                else
                    sc_objective = w_sigma * sigma_v + w_nu * norm(nu_v(:), 1);
                end
                
                minimize(sc_objective)
                
                cvx_end
                
                % Handle solver output
                if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
                    obj.X = X_v;
                    obj.U = U_v;
                    obj.nu = nu_v;
                    obj.sigma = sigma_v; % NEW: Save solved time
                    
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