function [ts, x, u] = sys_fwd(f, x0, ts, x_nom, u_hat, options, nvp)
%sys_fwd
%
% Forward integrator enforcing mixed constraints.
arguments
    f
    x0
    ts
    x_nom
    u_hat
    options
    nvp.initialization = false;
    nvp.beta1 = [];
    nvp.beta2 = [];
    nvp.db = [];
    nvp.g = [];
    nvp.gb = [];
    nvp.hessMinSolver = [];
    nvp.Vx = [];
    nvp.Vxx = [];
    nvp.Vxb = [];
end

x = zeros(length(ts), length(x0));
u = zeros(length(ts), length(u_hat(1,:)));
nx = length(x0);
nu = length(u_hat(1,:));
npsi = length(nvp.db);

% Initial conditions
x(1, :) = x0;

% Integration routine
for n = 1:length(ts)
    xx = x(n, :);
    xx = xx(:);
    xx_nom = x_nom(n,:)';
    uu_hat = u_hat(n,:)';

    % Get new controls
    u(n, :) = get_u;

    %% Advance the simulation
    if n<length(ts)
        % Timestep
        h = ts(n+1) - ts(n);
        switch options.integrator
            case "rk4"
                % Integrate one step fwd
                k_1 = full( f.fun( xx, u(n, :), ts(n) ) );
                k_2 = full( f.fun( xx + k_1*h/3, u(n, :), ts(n) + h/3 ) );
                k_3 = full( f.fun( xx - k_1*h/3 + k_2*h, u(n, :), ts(n) + h*2/3 ) );
                k_4 = full( f.fun( xx + k_1*h - k_2*h + k_3*h, u(n, :), ts(n) + h ) );
                x(n+1, :) = xx + (1/8)*( k_1 + 3*k_2 + 3*k_3+k_4 )*h;
            case "fwdEuler"
                k_1 = full( f.fun( xx, u(n, :), ts(n) ) );
                x(n+1, :) = xx + k_1*h;
        end
    end
    
end

    function u_new = get_u
        % Get new controls
        
        if nvp.initialization
            % State traj initialization: u_hat is u_nom
            u_new = uu_hat';
        elseif options.compTrick
            % Use the computational trick
            dx = ( xx - xx_nom );
            Vx = reshape(nvp.Vx(n,:), nx, 1);
            Vxx = reshape(nvp.Vxx(n,:), nx, nx);
            if isempty(nvp.Vxb)
                sol = nvp.hessMinSolver('x0', uu_hat, 'p', [xx; Vx + Vxx*dx; ts(n)], 'lbx', nvp.gb{3}(1), 'ubx', nvp.gb{3}(2), ...
                    'lbg', -inf, 'ubg', 0);
            else
                Vxb = reshape(nvp.Vxb(n,:), nx, length(nvp.db));
                sol = nvp.hessMinSolver('x0', uu_hat, 'p', [xx; Vx + Vxx*dx + Vxb*nvp.db; ts(n)], 'lbx', nvp.gb{3}(1), 'ubx', nvp.gb{3}(2), ...
                    'lbg', -inf, 'ubg', 0);
            end
            u_new = full(sol.x);
        else
            % Use beta1 and beta2
            if isempty(nvp.beta2)
                u_new = uu_hat' + ( reshape(nvp.beta1(n,:,:), nu, nx) * ( xx - xx_nom ) )';
            else
                u_new = uu_hat' + ( reshape(nvp.beta1(n,:,:), nu, nx) * ( xx - xx_nom ) )' + (reshape(nvp.beta2(n,:,:), nu, npsi) * nvp.db)';
            end
        end

        % Enforce constraints on u(t)
        % General constraints
        if ~isempty(nvp.g)
            % Check constraints violation
            constrValue = cellfun( @(c) full( c.fun(xx, u_new, ts(n)) ), nvp.g );
            activeConstraints = constrValue > 0;
            % If the constraints were violated, set u to satisfy them
            if any(activeConstraints)
                [~, activeConstraints] = max(constrValue);
                u_new =  fsolve(@(u) ...
                    full( nvp.g{activeConstraints}.fun(xx, u, ts(n)) ), ...
                    0, optimoptions('fsolve', 'Display', 'none') );
            end
        end
        % Box constraints
        if ~isempty(nvp.gb)
            % Enforce box constraints
            for m = 1:nu
                u_new(m) = min(max(u_new, nvp.gb{3}(m,1)), nvp.gb{3}(m,2));
            end
        end
    end


end