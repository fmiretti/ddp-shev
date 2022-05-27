% Series HEV
clear
close all
addpath(fullfile("..", "diffDP"))
% Import CasADi
addpath(fullfile(string(userpath), "Add-Ons", "casadi"))
import casadi.*

% Defined x,u,Vx,t as CasADi symbolic variables
x = MX.sym('x', 1);
Vx = MX.sym('Vx', 1);
u = MX.sym('u', 1);
t = MX.sym('t', 1);

%% Options
% Cosmetics
options.StateName = ["\sigma", "T_b"];
options.StateUnit = ["-", "°C"];
options.ControlName = "\tau_{gen}";
options.ControlUnit = "-";
options.CostUnit = "kg";
options.AdjointUnit = ["kg", "kg °C"];

% Max number of iterations
options.maxIter1 = 20;
options.maxIter2 = 10;
% Cost difference for convergence
options.eta1 = .001; % kg
options.eta2 = .001;
% Draw plots?
options.plots = true;
% Others
options.integrator = "rk4";

%% Define state dynamics, cost, time interval, nominal control trajectory and multipliers
[ts, x0, xf, f, L, F, u_bound, g, u_nom, b_nom] = shevProblemSetup("base");

%% Run
[x_nom, u_nom, V0, info] = diffDP(f, L, F, g, u_bound, x0, xf, ts, u_nom, b_nom, options);

% Evaluate the cumulative lagrangian
L0 = V0(end);

fprintf("Fuel consumption: %.1f g\n", L0*1e3)
fprintf("Final SOC: %.4f\n", x_nom(end,1))

