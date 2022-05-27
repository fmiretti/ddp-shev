function [time, x0, xf, f, L, F, u_bound, g, u_nom, b_nom] = shevProblemSetup(experiment, dynaprog)
%shevProblemSetup(experiment, dynaprog)
%
% Run the differential dynamic programming algorithm.
%
% Parameters
% ----------
% experiment : {'base', 'smallBatt', 'twoState'}
%   Experiment to run (See Ch. 7 of my PhD thesis).
% dynaprog : logical, default=False
%   Set to True for preparing a DynaProg experiment.
%
% Returns
% -------
% time : array
%   Time vector.
% x0 : array
%   Initial state.
% xf : array
%   Final state.
% f, L, F: function_handle
%   State dynamics, running cost and terminal cost.
% u_bound : array
%   Box constraints on the controls.
% g : function_handle cell array
%   Inequalty constraints in a cell array of function handles.
% u_nom : function_handle
%   Nominal control trajectory.
% b_nom : vector
%   Nominal endpoint multipliers.

arguments
    experiment {mustBeMember(experiment, {'base', 'smallBatt', 'twoState'})} = 'base' 
    dynaprog = false
end

%% Boundary conditions
% Environment temperature (°C)
switch experiment
    case 'twoState'
        Tenv = -10;
    otherwise
        Tenv = 20;
end

% Initial state
%   must be a column vector
switch experiment
    case 'twoState'
        x0 = [0.6; Tenv];
    otherwise
        x0 = 0.6;
end

% Final state
%   must be a column vector
xf = 0.6;

%% Load cycle data
load('WLTP3');
time = seconds(WLTP3.time); % s
vehSpd = WLTP3.vehSpd; % m/s
vehAcc = WLTP3.vehAcc; % m/s^2

%% Powertrain data
switch experiment
    case 'smallBatt'
        load("shev_data_smallBatt.mat")
    otherwise
        load("shev_data.mat")
end

%% Driveline model
for n = 1:length(time)
    % Wheels
    % Wheel speed (rad/s)
    wheelSpd  = vehSpd(n) ./ wh.radius;
    % Tractive Force (N)
    vehResForce =  veh.f0 + veh.f1 .* vehSpd(n) + veh.f2 .* vehSpd(n).^2;
    vehForce = (vehSpd(n)~=0) .* (vehResForce + veh.mass.*vehAcc(n));
    % Wheel torque (Nm)
    wheelTrq = vehForce .* wh.radius;

    motTrq(n) = ( wheelTrq < 0 ) .* mot.brakFactor .* wheelTrq ./ fd.spdRatio + ...
        ( wheelTrq >= 0 ) .* wheelTrq ./ fd.spdRatio;
    motSpd(n) = wheelSpd .* fd.spdRatio;

    motPwr(n) = mot.effMap( motSpd(n), motTrq(n) ) .* motSpd(n) .* motTrq(n);
end

if ~dynaprog
    % Use CasADi interpolant
    motPwr = casadi.interpolant('motPwr', 'linear', {time}, motPwr );
    % Convert MATLAB interpolant to CasADi interpolant
    batt.ocVolt = casadi.interpolant('ocVolt', 'linear', batt.ocVolt.GridVectors, batt.ocVolt.Values);
else
    % Use MATLAB interpolant
    motPwr = griddedInterpolant(time, motPwr);
end

%% State dynamics
%   must be a column vector function of x,u,t
SOC_lo = 0.4;
SOC_up = 0.8;

battPwr = @(u,t)  motPwr(t) - u .* gen.maxPwr;
switch experiment
    case {'base', 'smallBatt'}
        batt.eqRes = batt.eqRes(Tenv);
        battCurr = @(x,u,t) ( batt.ocVolt(x(1)) - sqrt( batt.ocVolt(x(1)).^2 - 4.*batt.eqRes.*battPwr(u,t) ) ) ./ ( 2*batt.eqRes );
        battCurr = memoize( battCurr );
        f_soc = @(x,u,t) - battCurr(x,u,t) ./ ( batt.nomCap * 3600 );
        if ~dynaprog
            % vanishing f above the upper soc threshold
            f = @(x,u,t) f_soc(x,u,t) .* ( 1 ./ ( 1 + exp(-4e2 * (0.8 - x)) ) ) .* ( 1 ./ ( 1 + exp(-4e2 * (x - 0.4)) ) ) ;
        else
            f = f_soc;
        end
    case 'twoState'
        batt.envTemp = 20;
        % Temperature in degrees C
        battCurr = @(x,u,t) ( batt.ocVolt(x(1)) - sqrt( batt.ocVolt(x(1)).^2 - 4 .* batt.eqRes(x(2)) .* battPwr(u,t) ) ) ./ ( 2*batt.eqRes(x(2)) );
        battCurr = memoize( battCurr );
        f_soc = @(x,u,t) - battCurr(x,u,t) ./ ( batt.nomCap * 3600 );
        f_temp = @(x,u,t) ( batt.heatTransCoef * ( batt.envTemp - x(2) ) + batt.eqRes( x(2) ) .* battCurr(x,u,t) .^ 2 ) ./ ( batt.heatCap );
        if ~dynaprog
            f = @(x,u,t) [f_soc(x,u,t); f_temp(x,u,t)];
        else
            battCurr = @(x,u,t) ( batt.ocVolt(x{1}) - sqrt( batt.ocVolt(x{1}).^2 - 4 .* batt.eqRes(x{2}) .* battPwr(u,t) ) ) ./ ( 2*batt.eqRes(x{2}) );
            f_soc = @(x,u,t) - battCurr(x,u,t) ./ ( batt.nomCap * 3600 );
            f_temp = @(x,u,t) ( batt.heatTransCoef * ( batt.envTemp - x{2} ) + batt.eqRes( x{2} ) .* battCurr(x,u,t) .^ 2 ) ./ ( batt.heatCap );
            f = {@(x,u,t) f_soc(x,u,t); @(x,u,t) f_temp(x,u,t)};
        end
end

%% Cost
% Running Cost
%   must be a scalar function of x,u,t
L = @(x,u,t) gen.fuelConsumption(u) * 1e-3; % gen.fuelConsumption is in g/s

% Terminal cost
%   must be a scalar function of x,t
F = @(x, t) 1e1*(x(1) - xf(1)).^2; % balancing term

%% Constraints 
% Control bounds u lower bound, u upper bound
u_bound = [0, 1];

% Mixed state-control inequality constraints g(x,u,t) <= 0
% SOC Window
a1 = 100;
a2 = a1;
g{1} = @(x,u,t) f_soc(x,u,t) + a1 * ( x(1) - SOC_up );
g{2} = @(x,u,t) - f_soc(x,u,t) + a2 * ( SOC_lo - x(1) );

%% Nominal multipliers
b_nom = - 1; 

%% Nominal control trajectory
switch experiment
    case 'smallBatt'
        u_nom = @(t) full( 0.05 .* ones(size(t)) .* (t<1545)  + 0.17 .* ones(size(t)) .* (t>=1545));
    otherwise
        u_nom = @(t) full( 0.06 .* ones(size(t)) );
end

end

