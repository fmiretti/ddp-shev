function defOptions = defaultOptions()
%defaultOptions() 
% 
% Get default options.
%  
% Returns
% -------
% struct
%   Options structure with default settings.

% Max number of iterations
defOptions.maxIter1 = 10;
defOptions.maxIter2 = 10;

% Convergence tolerances
defOptions.eta1 = .1;
defOptions.eta2 = .001;

% Enable computational trick
defOptions.compTrick = true;

% Neglect Hxx term in the base equations. May help preventing unbounded
%   solutions but generally slows down convergence.
defOptions.nullHxx = false;

% NLP solver options
defOptions.NLPmultiStart = false;

% Integrator
defOptions.integrator = "rk4";

% Analytical minimization. If empty: not available, use NLP solver.
defOptions.analytical_u_hat = [];

% Plot labels
defOptions.StateName = [];
defOptions.StateUnit = [];
defOptions.ControlName = [];
defOptions.ControlUnit = [];
defOptions.CostUnit = [];

% Draw plots?
defOptions.plots = true;
defOptions.plotTiling = "horizontal";

end