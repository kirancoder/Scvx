% load_parameters.m

% --- Global Parameters ---
K = 30;

% Max solver iterations for Successive Convexification
max_iterations = 50;

% Solver settings
solver_choice = 'SDPT3'; 
verbose_solver = false;

% --- Weight Constants ---
% Penalty weight for the virtual control
w_nu = 1e5; 

% --- Trust Region Parameters ---
% Initial trust region radius
tr_radius = 5.0;

% Trust region threshold variables for accepting/rejecting solutions
rho_0 = 0.0;
rho_1 = 0.25;
rho_2 = 0.9;

% Multipliers for shrinking (alpha) or expanding (beta) the trust region
alpha_tr = 2.0;
beta_tr = 3.2;

% Penalty weight for the total flight time
weight_sigma = 1.0; 
