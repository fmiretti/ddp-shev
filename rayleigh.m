% Example # 6.4  From Part 3, Sec. 6 of Martensson 1972

clear
close all

%% Options
% Max number of iterations
options.maxIter1 = 10;
options.maxIter2 = 10;
% Cost difference for convergence
options.eta1 = .1;
options.eta2 = .01;
% Draw plots?
options.plots = true;

%% Boundary conditions
% Initial state
%   must be a column vector
x0 = [-5; -5];
% Final state
%   must be a column vector
xf = [1; 3];

%% State dynamics
%   must be a column vector function of x,u,t
f = @(x,u,t) [x(2); -x(1) + 1.4*x(2) - 0.14*x(2)^3 + 4*u];

%% Cost
% Running Cost
%   must be a scalar function of x,u,t
L = @(x,u,t) x(1).^2 + u.^2;

% Terminal cost
%   must be a scalar function of x,t
F = @(x, t) 0;

%% Constraints 
% Control bounds u lower bound, u upper bound
u_bound = [-1, 1];

% Mixed state-control inequality constraints g(x,u,t) <= 0
g = {};

%% Time interval
t0 = 0;
tf = 2.5;
% Discretization timesteps for integrations
num_timesteps = 400;

ts = linspace(t0, tf, num_timesteps);

%% Inital values
% Inital lagrange multipliers
b = [0; 0]; 

% Initial nominal control trajectory
u_nom = @(t) interp1([t0  tf], [-0.5 -0.5], t);

%% Run
[x_nom, u_nom, V0, info] = diffDP(f, L, F, g, u_bound, x0, xf, ts, u_nom, b, options);
