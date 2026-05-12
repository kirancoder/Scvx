% main_2D_landing.m
% MATLAB implementation of the Free-Final-Time SCvx for 2D Rocket Landing.

clc;
clear;
close all;

% 1. Load parameters into the workspace (Make sure weight_sigma is defined here!)
load_parameters;

% 2. Initialize the model (2D Rocket Landing)
m = RocketLanding2D();
m.nondimensionalize();

% 3. Initialize trajectory using linear interpolation
[X, U] = m.initialize_trajectory(K);
sigma = m.t_f_guess;

% Arrays to store the history for logging/plotting
all_X = cell(1, max_iterations + 1);
all_U = cell(1, max_iterations + 1);
all_sigma = zeros(1, max_iterations + 1);

all_X{1} = m.x_redim(X);
all_U{1} = m.u_redim(U);
all_sigma(1) = sigma;

% 4. Setup Integrator and Problem wrapper (using the Free-Final-Time versions)
integrator = Discretization(m, K);
problem = SCProblem(m, K);

last_nonlinear_cost = [];
converged = false;

% 5. Start Successive Convexification Loop
for it = 1:max_iterations
    t0_it = tic;
    fprintf('--------------------------------------------------\n');
    fprintf('------------------ Iteration %02d ------------------\n', it);
    fprintf('--------------------------------------------------\n');

    t0_tm = tic;
    % NEW: calculate_discretization now requires sigma and returns S_bar
    [A_bar, B_bar, C_bar, S_bar, z_bar] = integrator.calculate_discretization(X, U, sigma);
    fprintf('Time for transition matrices: %.4f s\n', toc(t0_tm));

    % NEW: Pass S_bar, sigma_last, and weight_sigma to the CVX problem wrapper
    problem.set_parameters('A_bar', A_bar, 'B_bar', B_bar, 'C_bar', C_bar, 'S_bar', S_bar, 'z_bar', z_bar, ...
                           'X_last', X, 'U_last', U, 'sigma_last', sigma, ...
                           'weight_nu', w_nu, 'weight_sigma', weight_sigma, 'tr_radius', tr_radius, ...
                           'verbose_solver', verbose_solver);

    % Trust Region Inner Loop
    while true
        % Solve the Convex Subproblem
        error_flag = problem.solve();
        
        if error_flag
            fprintf('Solver failed (likely numerically ill-conditioned). Shrinking trust region...\n');
            tr_radius = tr_radius / alpha_tr;
            
            % Hard break to prevent infinite loops if completely stuck
            if tr_radius < 1e-5
                error('Trust region shrank too much, solver is completely stuck.');
            end
            
            problem.set_parameters('tr_radius', tr_radius);
            continue; % Retry the CVX solve with the smaller radius
        end

        % Retrieve solutions
        new_X = problem.get_variable('X');
        new_U = problem.get_variable('U');
        nu_var = problem.get_variable('nu');
        new_sigma = problem.get_variable('sigma'); % NEW: Retrieve optimized time

        % Integrate nonlinear dynamics to check linearization accuracy (using new sigma)
        X_nl = integrator.integrate_nonlinear_piecewise(new_X, new_U, new_sigma);

        % Calculate Costs (L1 Norms)
        linear_cost_dynamics = norm(nu_var(:), 1);
        nonlinear_cost_dynamics = norm(new_X(:) - X_nl(:), 1);

        linear_cost_constraints = m.get_linear_cost();
        nonlinear_cost_constraints = m.get_nonlinear_cost(new_X, new_U);

        linear_cost = linear_cost_dynamics + linear_cost_constraints; % J
        nonlinear_cost = nonlinear_cost_dynamics + nonlinear_cost_constraints; % L

        if isempty(last_nonlinear_cost)
            last_nonlinear_cost = nonlinear_cost;
            X = new_X;
            U = new_U;
            sigma = new_sigma;
            break;
        end

        actual_change = last_nonlinear_cost - nonlinear_cost;   % delta_L
        predicted_change = last_nonlinear_cost - linear_cost;   % delta_J

        fprintf('\nFlight Time (sigma):  %.4f s\n', new_sigma);
        fprintf('Virtual Control Cost: %.4e\n', linear_cost_dynamics);
        fprintf('Constraint Cost:      %.4e\n\n', linear_cost_constraints);
        fprintf('Actual change:        %.4e\n', actual_change);
        fprintf('Predicted change:     %.4e\n\n', predicted_change);

        % Check for convergence
        if abs(predicted_change) < 1e-4
            converged = true;
            break;
        else
            rho = actual_change / predicted_change;
            if rho < rho_0
                % Reject solution, shrink trust region
                tr_radius = tr_radius / alpha_tr;
                fprintf('Trust region too large. Solving again with radius=%.4f\n', tr_radius);
            else
                % Accept solution
                X = new_X;
                U = new_U;
                sigma = new_sigma;
                fprintf('Solution accepted.\n');

                if rho < rho_1
                    fprintf('Decreasing radius.\n');
                    tr_radius = tr_radius / alpha_tr;
                elseif rho >= rho_2
                    fprintf('Increasing radius.\n');
                    tr_radius = tr_radius * beta_tr;
                end

                last_nonlinear_cost = nonlinear_cost;
                break;
            end
        end
        
        % Update trust region for next inner loop attempt
        problem.set_parameters('tr_radius', tr_radius);
        fprintf('--------------------------------------------------\n');
    end

    fprintf('\nTime for iteration: %.4f s\n\n', toc(t0_it));

    % Save dimensionalized variables for plotting
    all_X{it + 1} = m.x_redim(X);
    all_U{it + 1} = m.u_redim(U);
    all_sigma(it + 1) = sigma;

    if converged
        fprintf('Converged after %d iterations.\n', it);
        % Trim unused pre-allocation
        all_X = all_X(1:it+1);
        all_U = all_U(1:it+1);
        all_sigma = all_sigma(1:it+1);
        break;
    end
end

if ~converged
    fprintf('Maximum number of iterations reached without convergence.\n');
end

% --- Visualization ---
fprintf('Plotting final trajectory...\n');
final_X = all_X{end};
final_U = all_U{end};
final_sigma = all_sigma(end);

figure('Name', 'Free-Final-Time Rocket Landing 2D', 'Position', [100, 100, 600, 800]);
plot(final_X(1, :), final_X(2, :), 'Color', [0.8 0.8 0.8], 'LineWidth', 2);
hold on;
grid on;
xlabel('X (East) [m]');
ylabel('Y (Up) [m]');
title(sprintf('2D Landing - Iteration %d - Final Time: %.2f s', length(all_X)-1, final_sigma));

% Add attitude and thrust vectors
attitude_scale = 0.5;
thrust_scale = 0.05;

for k = 1:K
    rx = final_X(1, k);
    ry = final_X(2, k);
    theta = final_X(5, k);
    gimbal = final_U(1, k);
    thrust = final_U(2, k);
    
    % Attitude vector
    dx = sin(theta) * attitude_scale;
    dy = cos(theta) * attitude_scale;
    quiver(rx, ry, dx, dy, 0, 'b', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
    
    % Thrust vector
    Fx = -sin(theta + gimbal) * thrust * thrust_scale;
    Fy = -cos(theta + gimbal) * thrust * thrust_scale;
    quiver(rx, ry, Fx, Fy, 0, 'r', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
end
axis equal;