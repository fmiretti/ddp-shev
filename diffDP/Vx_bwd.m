function [ts, y, u_hat, beta1, Z, A] = Vx_bwd(ts, yf, x_nom, u_nom, hessMinSolver, f, H, g, gb, nx, nu, options)
%Vx_bwd
% 
% Integrate base equations for a, Vx, Vxx.

% Storage space for u_star(t), beta1(t), A = (Huu + mu*guu_hat)^-1 and Z.
nRows = length(ts);
u_hat = zeros(nRows, nu);
beta1 = zeros(nRows, nu, nx);
A = zeros(nRows, nu, nu);
Z = zeros(nRows, nu, nu);

% Storage space for y = [a, Vx(:), Vxx(:)];
y = zeros(nRows, length(yf));

% Initial conditions
y(end, :) = yf;

for n = length(ts):-1:1
    yy = y(n, :);
    yy = yy(:);

    xx_nom =  x_nom(n, :);
    xx_nom = xx_nom(:);
    uu_nom = u_nom(n, :);
    Vx = reshape(yy(2:(1+nx)), [1 nx]);
    tt = ts(n);

    % Check whether the solution is unbounded. If it is, there is no point
    % in carrying out the integration further.
    if any(isnan(y(n, :))) 
        if n>1
            y(n-1, :) = nan(size(y(n, :)));     
        end
       u_hat(n,:) = nan(nu);
       beta1(n,:,:) = nan(nu,nx);
       Z(n,:,:) = nan(nu);
       A(n,:,:) = nan(nu);
       continue
    end

    % Get hamiltonian-minimizing controls
    [uu_hat, g_hat] = get_u_hat;
    u_hat(n,:) = uu_hat;
        
    [k_1, beta1(n,:,:), Z(n,:,:), A(n,:,:)] = ...
        Vx_odefun( ts(n), yy, xx_nom, uu_nom, uu_hat, f, H, g_hat, nx, nu, options);
    if n>1
        % Timestep
        h = ts(n-1) - ts(n);
        switch options.integrator
            case "rk4"
                % Fourth order runge kutta
                k_2 = ...
                    Vx_odefun( ts(n) + h/2, yy + k_1*h/2, xx_nom, uu_nom, uu_hat, f, H, g_hat, nx, nu, options);
                k_3 = ...
                    Vx_odefun( ts(n) + h/2, yy + k_2*h/2, xx_nom, uu_nom, uu_hat, f, H, g_hat, nx, nu, options);
                k_4 = ...
                    Vx_odefun( ts(n) + h, yy + k_3*h, xx_nom, uu_nom, uu_hat, f, H, g_hat, nx, nu, options);
                y(n-1, :) = yy + (1/6)*( k_1 + 2*k_2 + 2*k_3 + k_4 )*h;
            case "fwdEuler"
                y(n-1, :) = yy + k_1*h;
        end
    end

end

    function [uu_hat, g_hat] = get_u_hat
    % Get H-minimizing controls

        if isempty(options.analytical_u_hat)
            % Solve the NLP

            % Use a three-point multistart?
            if options.NLPmultiStart
                u_guess = [gb{3}(1) + (gb{3}(2)-gb{3}(1))*0.05, uu_nom, gb{3}(2) - (gb{3}(2)-gb{3}(1))*0.05];
            else
                u_guess = uu_nom;
            end
            sol = hessMinSolver('x0', u_guess, 'p', [xx_nom; Vx(:); tt], 'lbx', gb{3}(1), 'ubx', gb{3}(2), ...
                'lbg', -inf, 'ubg', 0);
            [~, idx] = min( full(sol.f) );
            uu_hat = full( sol.x(idx) );
        else
            % An analytical expression is available
            uu_hat = full( options.u_hat(xx_nom, Vx, tt) );
            uu_hat = min( max(uu_hat, gb{3}(1)), gb{3}(2) );
        end

        %% Detect active constraints
        % Notes:
        %  - each control variable can activate at most one constraint at a
        %   time
        %  - this portion currently only supports scalar controls and is
        %   not very robust in general

        constrValue = cellfun( @(c) full( c.fun(xx_nom, uu_hat, tt) ), g );
        lbConstrValue = gb{3}(1) - uu_hat;
        ubConstrValue = uu_hat - gb{3}(2);

        activeConstraints = full( constrValue > -1e-6 );
        activeLbConstraints = full( lbConstrValue > -1e-2 );
        activeUbConstraints = full( ubConstrValue > -1e-2 );

        g_hat = {};
        if any([activeConstraints(:); activeLbConstraints(:); activeUbConstraints(:)])
            [~, activeConstraints] = max([constrValue(:); lbConstrValue(:); ubConstrValue(:)]);
            if activeConstraints <= length(g)
                g_hat = g(activeConstraints);
            elseif activeConstraints == length(g)+1
                g_hat{end+1} = gb{1};
            elseif activeConstraints == length(g)+2
                g_hat{end+1} = gb{2};
            end
        end

        if activeLbConstraints
            g_hat = gb(1);
        elseif activeUbConstraints
            g_hat = gb(2);
        end

    end

end
