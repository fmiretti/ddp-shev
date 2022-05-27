function [x_nom, u_nom, V0, info] = diffDP(fexpr, L, Fexpr, gexprs, u_bound, x0, xf, ts, u_nom, b, options)
%diffDP(fexpr, L, Fexpr, gexprs, u_bound, x0, xf, ts, u_nom, b, options)
%
% Run the differential dynamic programming algorithm.
%
% Parameters
% ----------
% fexpr : function_handle array
%   The state dynamics in a column array of function handles.
% L, Fexpr : function_handle
%   Running cost and terminal cost.
% gexprs : function_handle cell array
%   Inequalty constraints in a cell array of function handles.
% u_bound : array
%   Box constraints on the controls (lb <= u <= ub) in an array (scalar
%   control) or matrix (nonscalar controls), specified as [lb, ub].
% x0 : array
%   Initial states in a column array.
% xf : array
%   Terminal states in a column array. If non-constrained, leave empty.
%   If partially constrained, only specify the constrained states so that
%   length(xf) = number of constrained states.
% ts : array
%   Time vector.
% u_nom : array
%   Nominal control trajectory, specified at ts.
% b : vector
%   Nominal endpoint multipliers. Specify the constrained states only, as
%   for xf.
% options : struct
%   Options structure.
%
% Returns
% -------
% x_nom : array
%   Optimal state trajectory sampled ad ts.
% u_nom : array
%   Optimal control trajectory sampled ad ts.
% V0 : scalar
%   Optimal cost.
% info: struct
%   Additional solution information.
%
% See Also
% --------
% defaultOptions
arguments
    fexpr function_handle
    L function_handle
    Fexpr function_handle
    gexprs cell
    u_bound double
    x0 double
    xf double
    ts double
    u_nom {mustBeA(u_nom, ["double", "function_handle"])}
    b double
    options struct
end

%% Initialize nominal trajectories
tf = ts(end);
ts = ts(:);
if isa(u_nom, "function_handle")
    u_nom = u_nom(ts);
end

nx = length(x0);
nu = length(u_nom(1,:));

%% Set options
defOptions = defaultOptions();
defFields = string(fieldnames(defOptions));
for n = 1:length(defFields)
    if ~isfield(options, defFields(n))
        options.(defFields(n)) = defOptions.(defFields(n));
    end
end

%% Set up Casadi Stuff
import casadi.*

x = MX.sym('x', nx);
Vx = MX.sym('Vx', nx);
u = MX.sym('u', nu);
t = MX.sym('t', 1);

% State equations
fexpr = fexpr(x,u,t);
f.fun = Function('f', {x, u, t}, {fexpr}, {'x','u','t'}, {'xdot'});
% Create the Jacobian expressions
f.fx = Function('fx', {x, u, t}, {jacobian(fexpr, x)}, {'x','u','t'}, {'fx0'});
f.fu = Function('fu', {x, u, t}, {jacobian(fexpr, u)}, {'x','u','t'}, {'fu0'});

% Hamiltonian
Hexpr = L(x,u,t) + Vx.' * fexpr;
H.expr = Hexpr;
H.fun = Function('H', {x, u, Vx, t}, {Hexpr}, {'x','u','Vx','t'}, {'H0'});
% Create the Jacobian expressions
H.Hx = Function('Hx', {x, u, Vx, t}, {jacobian(Hexpr, x)}, {'x','u','Vx','t'}, {'Hx0'});
H.Hu = Function('Hu', {x, u, Vx, t}, {jacobian(Hexpr, u)}, {'x','u','Vx','t'}, {'Hu0'});
% Create the Hessian expressions
Hexpr = hessian(Hexpr, [x; u]);
H.Hxx = Function('Hxx', {x, u, Vx, t}, {Hexpr(1:nx, 1:nx)}, {'x','u','Vx','t'}, {'Hxx0'});
H.Hux = Function('Hux', {x, u, Vx, t}, {Hexpr(nx+1:end, 1:nx)}, {'x','u','Vx','t'}, {'Hux0'});
H.Huu = Function('Huu', {x, u, Vx, t}, {Hexpr(nx+1:end, nx+1:end)}, {'x','u','Vx','t'}, {'Huu0'});

% Terminal cost
Fexpr = Fexpr(x,t);
F.fun = Function('F', {x, t}, {Fexpr}, {'x','t'}, {'F0'});
% Create the Jacobian expressions
F.Fx = Function('Fx', {x, t}, {jacobian(Fexpr, x)'}, {'x','t'}, {'Hx0'});
% Create the Hessian expressions
Fexpr = hessian(Fexpr, x);
F.Fxx = Function('Fxx', {x, t}, {Fexpr}, {'x','t'}, {'Hxx0'});

% Constraints
g = cell(1,size(gexprs,2));
for n = 1:length(gexprs)
    gexpr = gexprs{1,n}(x,u,t);
    g{n}.expr = gexpr;
    g{n}.fun = Function(['g' num2str(n)], {x, u, t}, {gexpr}, {'x','u','t'}, {['g' num2str(n) '0']});
    % Create the Jacobian expressions
    g{n}.gx = Function(['g' num2str(n) 'x'], {x, u, t}, {jacobian(gexpr, x)}, {'x','u','t'}, {['g' num2str(n) 'x0']});
    g{n}.gu = Function(['g' num2str(n) 'u'], {x, u, t}, {jacobian(gexpr, u)}, {'x','u','t'}, {['g' num2str(n) 'u0']});
    % Create the Hessian expressions
    gexpr = hessian(gexpr, [x; u]);
    g{n}.gxx = Function(['g' num2str(n) 'xx'], {x, u, t}, {gexpr(1:nx, 1:nx)}, {'x','u','t'}, {['g' num2str(n) 'xx0']});
    g{n}.gux = Function(['g' num2str(n) 'ux'], {x, u, t}, {gexpr(nx+1:end, 1:nx)}, {'x','u','t'}, {['g' num2str(n) 'ux0']});
    g{n}.guu = Function(['g' num2str(n) 'uu'], {x, u, t}, {gexpr(nx+1:end, nx+1:end)}, {'x','u','t'}, {['g' num2str(n) 'uu0']});
end
% gb: bound constraints
%  g{1} lower threshold constraint function and derivatives
%  g{2} upper threshold constraint function and derivatives
%  g{3} numeric array with lower and upper threshold

gb = cell(2,1);
gbexpr = [u_bound(1) - u; u - u_bound(2)];
for n = 1:length(gb)
    gexpr = gbexpr(n);
    gb{n}.expr = gexpr;
    gb{n}.fun = Function(['gb' num2str(n)], {x, u, t}, {gexpr}, {'x','u','t'}, {['gb' num2str(n) '0']});
    % Create the Jacobian expressions
    gb{n}.gx = Function(['gb' num2str(n) 'x'], {x, u, t}, {jacobian(gexpr, x)}, {'x','u','t'}, {['gb' num2str(n) 'x0']});
    gb{n}.gu = Function(['gb' num2str(n) 'u'], {x, u, t}, {jacobian(gexpr, u)}, {'x','u','t'}, {['gb' num2str(n) 'u0']});
    % Create the Hessian expressions
    gexpr = hessian(gexpr, [x; u]);
    gb{n}.gxx = Function(['gb' num2str(n) 'xx'], {x, u, t}, {gexpr(1:nx, 1:nx)}, {'x','u','t'}, {['gb' num2str(n) 'xx0']});
    gb{n}.gux = Function(['gb' num2str(n) 'ux'], {x, u, t}, {gexpr(nx+1:end, 1:nx)}, {'x','u','t'}, {['gb' num2str(n) 'ux0']});
    gb{n}.guu = Function(['gb' num2str(n) 'uu'], {x, u, t}, {gexpr(nx+1:end, nx+1:end)}, {'x','u','t'}, {['gb' num2str(n) 'uu0']});
end
gb{3} = u_bound;

% Bounds-enforcing controls
if size(gexprs,1) > 1
    for n = 1:length(g)
        if ~isempty(gexprs{2,n})
            g{n}.sol = gexprs{2,n};
        end
    end
end

% Terminal state constraints
if isempty(xf)
    psiexpr = 0;
    b = 0;
else
    str = "x(" + find(~isnan(xf)) + ") - " + num2str( xf(~isnan(xf)) );
    str = "@(x,t) vertcat(" + strjoin(str, ", ") + ")";
    psiexpr = str2func(str);
    psiexpr = psiexpr(x,t);
end
psi.fun = Function('psi', {x, t}, {psiexpr}, {'x','t'}, {'psi0'});
% Create the Jacobian expressions
psi.psix = Function('psix', {x, t}, {jacobian(psiexpr, x)}, {'x','t'}, {'psix0'});
psi.psiu = Function('psiu', {x, t}, {jacobian(psiexpr, u)}, {'x','t'}, {'psiu0'});
npsi = length(b);

% Set up NLP solver for u_hat
if isempty(g)
    prob = struct('f', H.fun(x,u,Vx,t), 'p',  vertcat(x, Vx, t),  'x', u);
else
    gexpr = MX.zeros(length(g),1);
    for n = 1:length(g)
        gexpr(n) = g{n}.fun(x,u,t);
    end
    prob = struct('f', H.fun(x,u,Vx,t), 'p',  vertcat(x, Vx, t),  'x', u, 'g', gexpr );
end
opts.ipopt.print_level = 0;
opts.print_time = false;
opts.show_eval_warnings = false;
hessMinSolver = nlpsol('solver', 'ipopt', prob, opts);

%% Nominal state traj
[ts, x_nom, u_nom] = sys_fwd(f, x0, ts, zeros(length(ts), nx), u_nom, options, 'g', g, 'gb', gb, 'initialization', true);

%% Plots Setup
if options.plots
    % Set up layout
    nrows = lcm(nx, nu);
    sv_span = nrows/nx;
    cv_span = nrows/nu;
    if options.plotTiling == "vertical"
        tl = tiledlayout(nrows, 2);
        cv_span = [cv_span 1];
        sv_span = [sv_span 1];
    elseif options.plotTiling == "horizontal"
        tl = tiledlayout(2, nrows);
        cv_span = [1 cv_span];
        sv_span = [1 sv_span];
    end

    % Set default axes labels if unspecified
    if isempty(options.StateName)
        options.StateName = "x_" + (1:nx);
    end
    if isempty(options.StateUnit)
        options.StateUnit = strings(1,nx);
    end
    if isempty(options.ControlName)
        options.ControlName = "u_" + (1:nu);
    end
    if isempty(options.ControlUnit)
        options.ControlUnit = strings(1,nu);
    end

    % Create axes and axes labels
    for n = 1:nx
        ax(n) = nexttile(sv_span); %#ok<*AGROW> 
        if isempty(options.StateUnit(n)) || strcmp(options.StateUnit(n), "")
            ylabel(options.StateName(n))
        else
            ylabel(options.StateName(n) + ", " + options.StateUnit(n))
        end
    end
    for n = 1:nu
        ax(n+nx) = nexttile(cv_span); 
        if isempty(options.ControlUnit(n)) || strcmp(options.ControlUnit(n), "")
            ylabel(options.ControlName(n))
        else
            ylabel(options.ControlName(n) + ", " + options.ControlUnit(n))
        end
    end
    legend('Location', 'southwest')

    % Style the layout and the axes
    tl.TileSpacing = 'compact';
    tl.Padding = 'compact';
    for n = 1:length(ax)
        hold(ax(n), 'on')
        grid(ax(n), 'on')
        ax(n).XLabel.String = 'Time, s';
    end
    linkaxes(ax, 'x')
    xlim([ts(1), ts(end)])

    % Line styles
    lineColors = {'#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',...
        '#8c564b', '#7f7f7f', '#bcbd22', '#17becf'};
    linestyles = ["-", repmat(["--","-.","."], 1, options.maxIter1)];
end

%% Main loop
info.b_iter(1,:) = b;
info.endpointError = [];
info.F_iter = [];
info.Vx = {};

for iter2 = 1:options.maxIter2
    for iter1 = 1:options.maxIter1
        %% Live plots
        if options.plots
            if iter1 == 1
                linewidth = 1.5;
                mkrSize = 2;
            else
                linewidth = 0.5;
                mkrSize = 1;
            end
            for n = 1:nx
                plot(ax(n), ts, x_nom(:,n), linestyles(iter1), ...
                    'Color', lineColors{mod(iter2-1,9)+1}, 'LineWidth', linewidth, ...
                    'MarkerSize', mkrSize, ...
                    'DisplayName', sprintf("iter #%d.%d", iter2-1, iter1-1))
            end
            for n = 1:nu
                plot(ax(n+nx), ts, u_nom(:,n), linestyles(iter1), ...
                    'Color', lineColors{mod(iter2-1,9)+1}, 'LineWidth', linewidth, ...
                    'MarkerSize', mkrSize, ...
                    'DisplayName', sprintf("iter #%d.%d", iter2-1, iter1-1))
            end
            drawnow
        end

        %% Diagnostics
        % Verbose diagnostics
        Ls = zeros(length(ts), 1);
        for n = 1:length(ts)
            Ls(n) = full( L( x_nom(n,:), u_nom(n,:), ts(n) ) );
        end
        fprintf("iter #%d.%d\n", iter2-1, iter1-1)
        fprintf("Cumulative lagrangian: %.2f;\t Terminal cost: %.2f;\t FXEP error: %.2f.\n", ...
            trapz(ts, Ls), full( F.fun( x_nom(end,:), ts(end) ) ), b.' * full( psi.fun(x_nom(end,:), ts(end)) ))

        % Save new cost
        info.V0_aug_iter(iter1,iter2) = totalCost(x_nom, u_nom, ts, L, F, b, psi);
        info.V0_iter(iter1,iter2) = totalCost(x_nom, u_nom, ts, L, struct('fun', @(x,t) 0 ), 0, psi);
        if ~isempty(xf)
            info.endpointError(end+1,:) = ( x_nom(end,1:npsi) - xf' );
        end
        info.F_iter(end+1,:) = full( F.fun(x_nom(end,:), tf));

        %% Integrate base equations
        % Evaluate boundary conditions
        a_f = 0;
        Vx_f = full( F.Fx( x_nom(end,:), ts(end) ) + psi.psix( x_nom(end,:), ts(end) ) .' *  b );
        Vxx_f = full( F.Fxx( x_nom(end,:), ts(end) ) ) + tensor_times_vector(zeros(npsi,nx,nx), b, 1);
        yf = [a_f; Vx_f(:); Vxx_f(:)];

        % Integrate a, Vx, Vxx backwards in time
        % ys: a, Vx, Vxx
        [ts, ys, u_hat, beta1, Z, A] = Vx_bwd(ts, yf, x_nom, u_nom, ...
            hessMinSolver, f, H, g, gb, nx, nu, options);

        a = ys(:,1);
        Vx = ys(:,2:nx+1);
        Vxx = ys(:,nx+2:end);
        
        % Store Vx (adjoint variables) 
        info.Vx{end+1} = {Vx};

        % Time t_eff
        N_eff = find(abs(a) < options.eta1, 1, 'first');
        
        % Check convergence
        if N_eff <= 1
            % Break
            fprintf("Solved free endpoint problem at iter #%d\n\n", iter1)
            break
        end
        if iter1 == options.maxIter1
            fprintf("Didn't solve free endpoint problem because the " + ...
                "maximum number of iterations (%d) was reached.\n\n", iter1);
            break
        end

        % Jacobson Step-Size Adjustment
        if options.compTrick
            [x_nom, u_nom, halt] = stepSizeAdj(f, L, F, psi, ys, ts,...
                x_nom, u_nom, u_hat, b, N_eff, options, ...
                "g", g, "gb", gb, "hessMinSolver", hessMinSolver);
        else
            [x_nom, u_nom, halt] = stepSizeAdj(f, L, F, psi, ys, ts, ...
                x_nom, u_nom, u_hat, b, N_eff, options, ...
                "beta1", beta1, "g", g, "gb", gb);
        end
        
        if halt
            break
        end

    end

    %% Check convergence (fixed endpoint problem)
    if full( norm(psi.fun(x_nom(end,:), tf)) ) < options.eta2 || options.maxIter2 == 1
        % Break
        fprintf("Converged at iter #%d\n\n", iter2-1)
        break
    end
    if iter2 == (options.maxIter2+1)
        fprintf("Ended at iter #%d because the maximum number of iterations was reached.\n\n", iter)
        break
    end

    %% Integrate base equations
    % Evaluate boundary conditions
    Vb_f = full( psi.fun(x_nom(end,:), tf).' );
    Vxb_f = full( psi.psix(x_nom(end,:), tf) );
    Vbb_f = zeros(npsi);
    yf = [Vb_f(:); Vxb_f(:); Vbb_f(:)];

    % Integrate Vxb,Vbb backwards in time
    % ys: Vb, Vxb, Vbb
    [~, ys, beta2] = Vb_bwd(ts, yf, x_nom, u_nom, u_hat, beta1, A, Z, f, nx, nu, npsi, options);

    % Multipliers adjustment method
    if options.compTrick
        [x_nom, u_nom, b] = mltpAdj(f, L, F, psi, ys, ts, ...
            x_nom, u_nom, u_hat, b, options, ...
            "a", a, "g", g, "gb", gb, ...
            "hessMinSolver", hessMinSolver, "Vx", Vx, "Vxx", Vxx);
    else
        [x_nom, u_nom, b] = mltpAdj(f, L, F, psi, ys, ts, ...
            x_nom, u_nom, u_hat, b, options, ...
            "a", a, "g", g, "gb", gb, "beta1", beta1, "beta2", beta2);
    end

    % Save new multipliers
    info.b_iter(iter2+1,:) = b;

end

%% Plot cost vs iterations
if options.plots
    costPlot(info.V0_iter, info.V0_aug_iter, info.endpointError, options)
end

%% Return info
V0 = info.V0_iter( find(info.V0_iter(:)~=0, 1, 'last') );

info.a = a;
info.Vxx = Vxx;
info.ts = ts;
info.options = options;

end
