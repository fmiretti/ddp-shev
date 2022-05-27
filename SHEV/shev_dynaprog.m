% Run experiments with dynamic programming 
% To run this script, you need to install DynaProg:
% https://www.mathworks.com/matlabcentral/fileexchange/84260-dynaprog

clear
%% Select experiment
experiment = "base";

%% Set up the problem
[ts, x0, xf, f, L, F, u_bound] = shevProblemSetup(experiment, true);

% Initial state 
x_init = num2cell(x0);

% State variable grid
switch experiment
    case 'twoState'
        SVnames = {'SOC', 'Battery Temperature, °C'};
        x_grid = {0.4:0.0005:0.8, -10:0.1:20};
        % Final state constraints
        x_final = {[xf xf+0.001], []};
    otherwise
        SVnames = 'SOC';
        x_grid = {0.4:0.0005:0.8};
        % Final state constraints
        x_final = {[xf xf+0.001]};
end

% Control variable grid
CVnames = {'Genset power'};
u1_grid = linspace(0, 1, 81);
u_grid = {u1_grid};

% Create exogenous input
w{1} = ts(1:end-1);
w{2} = ts(2:end) - ts(1:end-1);

% Number of stages (time intervals)
Nint = length(w{1});

fun = @(x, u, w) dp_model(x, u, w, f, L);

% Create DynaProg object
prob = DynaProg(x_grid, x_init, x_final, u_grid, Nint, ...
 fun, 'ExogenousInput', w);

%% Solve and visualize results
% Solve the problem
prob = run(prob);

%% Rerun the simulation with the rk4 integrator
u_dp_sol = prob.ControlProfile{1};
% Reload system dynamics
[~, ~, ~, f, L] = shevProblemSetup(experiment, false);
ff.fun = f;
nx = 2;
[ts, x_dp_sol, u_dp_sol] = sys_fwd(ff, x0, ts, zeros(length(ts), nx), [u_dp_sol(:); 0], struct("integrator", "rk4", "compTrick", false), "initialization", true);
cost_dp_sol = totalCost(x_dp_sol, u_dp_sol, ts, L, struct('fun', @(x,t) 0), 0, struct('fun', @(x,t) 0));

%% Plot
% Add SV, CV and cost names to be used in the plot
prob.StateName = SVnames;
prob.ControlName = CVnames;
prob.CostName = 'Fuel Consumption, g';
prob.Time = ts;

figure
t = plot(prob);

fprintf("Fuel consumption: %.1f g\n", cost_dp_sol*1e3)
fprintf("Final SOC: %.4f\n", x_dp_sol(end,1))

