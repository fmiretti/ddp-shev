clear
%% EM data
% Main data
mot.maxPwr = 85e3; % W
mot.maxSpd = 8000; % rpm
mot.maxSpd = mot.maxSpd * pi/30;
maxTrq = 300; % Nm

mot.inertia = 0;  % g*m^2
mot.mass = 10; % kg

% Efficiency map
trqBrkNorm = linspace(-1,1,11*2-1);
spdBrkNorm = linspace(0,1,14);
eff = [52 78 80 82 82 80 78 78 76 74 70;52 78 80 82 82 80 78 78 76 74 70;52 78 84 86 86 86 86 84 84 82 82;52 80 86 86 88 88 88 86 86 84 84;52 80 86 88 88 88 88 88 86 86 86;52 80 86 88 90 90 90 90 88 88 88;52 82 88 90 91 91 92 92 92 91 91;52 82 88 91 91 92 93 93 93 93 93;52 84 90 92 92 93 93 93 93 93 93;62 86 91 93 93 93 93 93 93 93 93;66 86 91 93 93 93 93 93 93 93 93;70 86 92 93 93 93 93 93 93 93 93;70 86 92 93 93 93 93 93 93 93 93;72 86 91 93 93 93 93 93 93 93 93];
eff = [eff(:,end:-1:2) eff];
eff = eff./100;

emTrqBrk = trqBrkNorm .* maxTrq;
emSpdBrk = spdBrkNorm .* mot.maxSpd;
mot.effMap = griddedInterpolant({emSpdBrk, emTrqBrk}, eff);

% Limit torque
kneeSpd = mot.maxPwr/maxTrq;

maxTrqData = maxTrq .* (emSpdBrk<kneeSpd) + mot.maxPwr ./ emSpdBrk .* (emSpdBrk>=kneeSpd);
maxTrqData(1) = maxTrq;

mot.maxTrq = griddedInterpolant(emSpdBrk, maxTrqData);
mot.minTrq = griddedInterpolant(emSpdBrk, -maxTrqData);

% Original data
% maxTrq = 200; % Nm
% maxPwr = 30e3; % W
% maxSpd = 6000; % rpm

%% Generator
eng = load("engData.mat");
eng = scaleEngineMap(eng, 61e3);

em = mot;
em.tcSpdRatio = 1.2;

% Number of discretized gen-set power values
numPoints = 40;
% Discretized gen-set power values
elPwr = linspace(0, eng.maxPwr, 40);
% Initial guesses for the engine speed
engSpd0 = linspace(eng.idleSpd, eng.maxSpd, 40);
% Lower and upper bounds for the engine speed and torque
lb = [eng.idleSpd; 0];
ub = [eng.maxSpd; inf];
% Solver options
options = optimoptions("fmincon",...
    "Algorithm","interior-point",...
    "EnableFeasibilityMode",true,...
    "SubproblemAlgorithm","cg");

% Run the solver
for n = 1:length(elPwr)
    [optPoint(n,:), ~, flag(n)] = fmincon(@(x) eng.fuelMap(x(1), x(2)), [engSpd0(n), 10], [], [], [], [], lb, ub, ...
        @(x) elPwrCon(x, eng, em, elPwr(n)), options);
end

% Remove failed points
optPoint(flag<1,:) = [];
elPwr(flag<1) = [];

%%
gen.maxPwr = max(elPwr);
pp = polyfit(elPwr ./ gen.maxPwr, eng.fuelMap(optPoint(:,1), optPoint(:,2)), 2);
gen.fuelCoefs = pp;
% gen.fuelConsumption = @(x) polyval(gen.fuelCoefs, x);
gen.fuelConsumption = @(x) pp(1) .* x.^2 + pp(2) .* x + pp(3);
gen.info = "fuelConsumption is in (g/s) and is expressed as a function " + ...
    "of normalized electric power, i.e. elPwr/gen.maxPwr";

%% Plots
close all
[fig, contourAx] = engMapPlot(eng);
colormap(contourAx, pyColorMap("viridis"))
scatterAx = axes(fig);
scatter(scatterAx, optPoint(:,1)*30/pi, optPoint(:,2), 50, elPwr ./ gen.maxPwr, ...
    'filled');
colormap(scatterAx, pyColorMap("plasma"))
% colormap(scatterAx, pyColorMap("inferno"))
% colormap(scatterAx, pyColorMap("magma"))
% colormap(scatterAx, pyColorMap("cividis"))

scatterAx.UserData = linkprop([contourAx, scatterAx],...
    {'Position','InnerPosition','xtick','ytick', ...
    'ydir','xdir','xlim','ylim'});
scatterAx.Visible = 'off';
lgd = findobj(fig.Children, 'Type', 'Legend');
lgd.Location = 'northwest';

exportFigureThesis(fig, 'Height', 11, 'FileName', 'engScatter')

%%
close all
[fig, contourAx] = emMapPlot(em);
colormap(contourAx, pyColorMap("viridis"))
scatterAx = axes(fig);
scatter(scatterAx, optPoint(:,1)*30/pi*em.tcSpdRatio, optPoint(:,2)/em.tcSpdRatio, 50, elPwr ./ gen.maxPwr, ...
    'filled');
colormap(scatterAx, pyColorMap("plasma"))

scatterAx.UserData = linkprop([contourAx, scatterAx],...
    {'Position','InnerPosition','xtick','ytick', ...
    'ydir','xdir','xlim','ylim'});
scatterAx.Visible = 'off';
lgd = findobj(fig.Children, 'Type', 'Legend');
lgd.Location = 'northwest';

exportFigureThesis(fig, 'Height', 10, 'FileName', 'emScatter')

%%
close all
figure
scatter(elPwr ./ gen.maxPwr, eng.fuelMap(optPoint(:,1), optPoint(:,2)));
hold on
fplot(@(x) polyval(pp, x), [-0.2 1.2]);

%% Final drive
fd.spdRatio = 4.2; % final drive
% fd.inertia = 0; % kg*m^2, rotational inertia of final drive, measured at input

%% Battery
% LiFePO
% main data based on
% thermal data based on 10.1016/j.jpowsour.2015.06.129

% Charcteristics
socBrk = 0.1:0.05:0.9;
cellOCV = [3.23740854764228;3.26575089639468;3.28674170326015;3.30599868708030;3.32504787988464;3.33346144837672;3.33871139954655;3.34478699733306;3.34938249667319;3.35431779549635;3.35957734516711;3.36684827133382;3.37553750385854;3.38517568202627;3.39257397772322;3.40418319672909;3.42434317055734];
cellReq = [0.00312541995262335;0.00295569040810214;0.00286378934122263;0.00283616676130315;0.00284618603238887;0.00287277790410816;0.00291074346415429;0.00295488380022049;0.00300000000000000;0.00304285559469203;0.00308919349833080;0.00314728387286201;0.00322539688023142;0.00333180268238473;0.00360476133631698;0.00409569001243483;0.00468259454280150];

% Other data
cellCap = 5.4; % Ah
batt.coulombic_eff = .99;  % -;

% Battery arrangement
Ns = 51;
Np = 2;

batt.eqRes = griddedInterpolant(socBrk, cellReq*Ns/Np);
batt.ocVolt = griddedInterpolant(socBrk, cellOCV*Ns);
batt.nomCap = cellCap*Np; % Ah
batt.nomEnergy = batt.nomCap * 3.4 * Ns; % Wh

% temp-dependent equivalent resistance
batt.eqRes = batt.eqRes(0.6);
batt.eqRes = @(T) batt.eqRes .* (0.0010 * T.^2 - 0.0903 .* T + 2.4520);

% Limit current
% minVolt = 2.6;
% maxVolt = 3.6;
% 
% minCurrData = ( batt.ocVolt.Values - maxVolt*Ns ) ./ batt.eqRes.Values;
% maxCurrData = ( batt.ocVolt.Values - minVolt*Ns ) ./ batt.eqRes.Values;
% 
% maxCRate = 18;
% minCurrData = max(minCurrData, - batt.nomCap*maxCRate);
% maxCurrData = min(maxCurrData, batt.nomCap*maxCRate);
% 
% batt.minCurr = griddedInterpolant(socBrk, minCurrData);
% batt.maxCurr = griddedInterpolant(socBrk, maxCurrData);

maxCRate = 10;
batt.minCurr = - batt.nomCap*maxCRate;
batt.maxCurr = batt.nomCap*maxCRate;

% Mass. Assume 20% additional mass for BMS and container.
energyDensity = 100; % Wh/kg
batt.mass = (batt.nomEnergy / energyDensity) * 1.2;

% Thermal data based on 10.1016/j.jpowsour.2015.06.129
%     Temperature in degrees C valid from 0-40 °C
%     batt.eqRes = @(T) batt.eqRes(0.6) .* (0.0010 * T.^2 - 0.0903 .* T + 2.4520);
batt.specHeatCap = 987; % J/kgK
batt.heatCap = batt.specHeatCap .* batt.mass ./ 1.2;
batt.heatTransCoef = 4.343; % W/K 

%% Vehicle Data
veh.mass = 1200; % kg m
veh.f0 = 76.11; % N
veh.f1 = 2.957 ; %  N/(m/s)
veh.f2 = 0.3664; % N/((m/s)^2)

veh.gravity = 9.81;    % m/s^2

%% Wheels
wh.radius = 0.3;    % (m)

%% Store data
save('noteData.mat', 'veh', 'wh', 'fd', 'gen', 'mot', 'batt')


function [c,ceq] = elPwrCon(x, eng, em,  elPwr)
    % x(1): engSpd, x(2): engTrq
    % Maximum torque (inequality) constraint
    c = x(2) - eng.maxTrq(x(1));
    % Power balance (equality) constraint
    ceq = elPwr - em.effMap(x(1) .* em.tcSpdRatio, x(2) ./ em.tcSpdRatio) .* x(1) .* x(2);
end