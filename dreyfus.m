% Example # 6.4  From Part 3, Sec. 6 of Martensson 1972

clear
close all

%% Options
% Max number of iterations
options.maxIter1 = 10;
options.maxIter2 = 20;
% Cost difference for convergence
options.eta1 = .01;
options.eta2 = .05;
% Draw plots?
options.plots = true;

%% Boundary conditions
% Initial state
%   must be a column vector
x0 = [0; 0];
% Final state
%   must be a column vector
xf = [1e5; 0];

%% State dynamics
%   must be a column vector function of x,u,t
a = 64; % ft/s2
g = 32; % ft/s
f = @(x,u,t) [x(2); a .* sin(u) - g];

%% Cost
% Running Cost
%   must be a scalar function of x,u,t
L = @(x,u,t) -a .* cos(u);

% Terminal cost
%   must be a scalar function of x,t
F = @(x, t) 0;

%% Constraints 
% Control bounds u lower bound, u upper bound
u_bound = [-4 4];

% Mixed state-control inequality constraints g(x,u,t) <= 0
g = {};

%% Time interval
t0 = 0;
tf = 100;
% Discretization timesteps for integrations
num_timesteps = 200;

ts = linspace(t0, tf, num_timesteps);

%% Inital values
% Inital lagrange multipliers
b = [0; 0]; 

% Initial nominal control trajectory
u_nom = @(t) 1.6 - 1.5 .* (t/100);

%% Run
[x_nom, u_nom, V0, info] = diffDP(f, L, F, g, u_bound, x0, xf, ts, u_nom, b, options);
