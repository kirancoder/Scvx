% main_free_3D_landing.m
% MATLAB implementation of 6-DoF 3D Rocket Landing via Free-Final-Time SCvx

clc;
clear;
close all;

% 1. Load parameters (Must include weight_sigma!)
load_parameters;

% 2. Initialize 3D Model
m = RocketLanding3D(K);
m.nondimensionalize();

% 3. Initialize trajectory
[X, U] = m.initialize_trajectory(K);
sigma = m.t_f_guess;

% Arrays to store the history
all_X = cell(1, max_iterations + 1);
all_U = cell(1, max_iterations + 1);
all_sigma = zeros(1, max_iterations + 1);

all_X{1} = m.x_redim(X);
all_U{1} = m.u_redim(U);
all_sigma(1) = sigma;

% 4. Setup Integrator and Problem wrapper (Free-Final-Time versions)
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
    % NEW: calculate_discretization handles sigma and S_bar
    [A_bar, B_bar, C_bar, S_bar, z_bar] = integrator.calculate_discretization(X, U, sigma);
    fprintf('Time for transition matrices: %.4f s\n', toc(t0_tm));

    % NEW: Pass S_bar, sigma_last, and weight_sigma to CVX
    problem.set_parameters('A_bar', A_bar, 'B_bar', B_bar, 'C_bar', C_bar, 'S_bar', S_bar, 'z_bar', z_bar, ...
                           'X_last', X, 'U_last', U, 'sigma_last', sigma, ...
                           'weight_nu', w_nu, 'weight_sigma', weight_sigma, 'tr_radius', tr_radius, ...
                           'verbose_solver', verbose_solver);

    while true
        error_flag = problem.solve();

        if error_flag
            fprintf('Solver failed (likely numerically ill-conditioned). Shrinking trust region...\n');
            tr_radius = tr_radius / alpha_tr;
            
            if tr_radius < 1e-4
                if ~isempty(last_nonlinear_cost)
                    fprintf('\nTrust region reached noise floor (%.5f). Forcing convergence!\n', tr_radius);
                    converged = true;
                    break; % Break the inner while loop
                else
                    error('Trust region shrank completely on the first iteration! Check initial guess.');
                end
            end
            
            problem.set_parameters('tr_radius', tr_radius);
            continue; 
        end

        new_X = problem.get_variable('X');
        new_U = problem.get_variable('U');
        nu_var = problem.get_variable('nu');
        new_sigma = problem.get_variable('sigma'); % Retrieve optimized time

        % Integrate nonlinear dynamics with new_sigma
        X_nl = integrator.integrate_nonlinear_piecewise(new_X, new_U, new_sigma);

        % Calculate Costs
        linear_cost_dynamics = norm(nu_var(:), 1);
        nonlinear_cost_dynamics = norm(new_X(:) - X_nl(:), 1);

        linear_cost_constraints = problem.model_linear_cost;
        nonlinear_cost_constraints = m.get_nonlinear_cost(new_X, new_U);

        linear_cost = linear_cost_dynamics + linear_cost_constraints;
        nonlinear_cost = nonlinear_cost_dynamics + nonlinear_cost_constraints;

        if isempty(last_nonlinear_cost)
            last_nonlinear_cost = nonlinear_cost;
            X = new_X;
            U = new_U;
            sigma = new_sigma;
            break;
        end

        actual_change = last_nonlinear_cost - nonlinear_cost;
        predicted_change = last_nonlinear_cost - linear_cost;

        fprintf('\nFlight Time (sigma):  %.4f s\n', new_sigma);
        fprintf('Virtual Control Cost: %.4e\n', linear_cost_dynamics);
        fprintf('Constraint Cost:      %.4e\n\n', linear_cost_constraints);
        fprintf('Actual change:        %.4e\n', actual_change);
        fprintf('Predicted change:     %.4e\n\n', predicted_change);

        if abs(predicted_change) < 1e-4
            converged = true;
            break;
        else
            rho = actual_change / predicted_change;
            if rho < rho_0
                tr_radius = tr_radius / alpha_tr;
                fprintf('Trust region too large. Solving again with radius=%.4f\n', tr_radius);
            else
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
        problem.set_parameters('tr_radius', tr_radius);
    end

    fprintf('\nTime for iteration: %.4f s\n\n', toc(t0_it));

    all_X{it + 1} = m.x_redim(X);
    all_U{it + 1} = m.u_redim(U);
    all_sigma(it + 1) = sigma;

    if converged
        fprintf('Converged after %d iterations.\n', it);
        all_X = all_X(1:it+1);
        all_U = all_U(1:it+1);
        all_sigma = all_sigma(1:it+1);
        break;
    end
end

if ~converged
    fprintf('Max iterations reached without convergence.\n');
end

% --- 3D Visualization ---
fprintf('Plotting final 3D trajectory...\n');
final_X = all_X{end};
final_U = all_U{end};
final_sigma = all_sigma(end);

figure('Name', 'Free-Final-Time Rocket Landing 3D', 'Position', [100, 100, 800, 800]);
plot3(final_X(2, :), final_X(3, :), final_X(4, :), 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
hold on; grid on; view(3);
xlabel('X (East) [m]'); ylabel('Y (North) [m]'); zlabel('Z (Up) [m]');
title(sprintf('3D Landing - Iteration %d - Final Time: %.2f s', length(all_X)-1, final_sigma));

% Draw target pad
th = linspace(0, 2*pi, 50);
fill3(20*cos(th), 20*sin(th), zeros(1,50), [0.9 0.9 0.9], 'EdgeColor', 'k');

attitude_scale = 30;
thrust_scale = 0.0001;

for k = 1:K
    rx = final_X(2, k);
    ry = final_X(3, k);
    rz = final_X(4, k);
    
    q = final_X(8:11, k);
    CBI = m.dir_cosine(q);
    
    % Direction mapping logic
    dx = CBI' * [0; 0; 1];
    Fx = CBI' * final_U(:, k);
    
    % Attitude vector
    quiver3(rx, ry, rz, dx(1)*attitude_scale, dx(2)*attitude_scale, dx(3)*attitude_scale, ...
            0, 'b', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
    
    % Thrust vector
    quiver3(rx, ry, rz, -Fx(1)*thrust_scale, -Fx(2)*thrust_scale, -Fx(3)*thrust_scale, ...
            0, 'r', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
end
axis equal;