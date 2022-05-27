% Example # 6.4  From Part 3, Sec. 6 of Martensson 1972

clear
close all
addpath(fullfile(string(userpath), "Add-Ons", "casadi"))

%% Options
% Max number of iterations
options.maxIter1 = 10;
options.maxIter2 = 10;
% Cost difference for convergence
options.eta1 = .002;
options.eta2 = .002;
% Draw plots?
options.plots = true;

%% Boundary conditions
% Initial state
%   must be a column vector
x0 = [0; 1];
% Final state
%   must be a column vector
xf = [0; -1];

%% State dynamics
%   must be a column vector function of x,u,t
mass = 1; % kg
f = @(x,u,t) [x(2); u(1)/mass];

%% Cost
% Running Cost
%   must be a scalar function of x,u,t
L = @(x,u,t) 0.5 * u(1).^2;

% Terminal cost
%   must be a scalar function of x,t
F = @(x,u,t) 20 * (x(:) - xf).' * (x(:) - xf);

%% Constraints 
% Control bounds u lower bound, u upper bound
u_bound = [-10, 10];

% Mixed state-control inequality constraints g(x,u,t) <= 0
a1 = 45;
a2 = 500;
g{1} = @(x,u,t) u + a1*x(2) + a2*(x(1) - 1/9);

%% Time interval
t0 = 0;
tf = 1;
% Discretization timesteps for integrations
num_timesteps = 200;

ts = linspace(t0, tf, num_timesteps);

%% Inital values
% Inital lagrange multipliers
b = [-5; 5];

% Initial nominal control trajectory
u_nom = @(t) interp1([t0  tf], [-5 -5], t);

%% Run
[x_nom, u_nom, V0, info] = diffDP(f, L, F, g, u_bound, x0, xf, ts, u_nom, b, options);

% Evaluate the cumulative lagrangian
L0 = V0(end);

fprintf("Cumulative Lagrangian: %.3f\n", L0)
fprintf("Final state: %.4f, %.4f\n", x_nom(end,:))
