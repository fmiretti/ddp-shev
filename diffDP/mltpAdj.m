function [x_nom, u_nom, b] = mltpAdj(f, L, F, psi, ys, ts, x_nom, u_nom, u_hat, b, options, nvp)
%mltpAdj(f, L, F, psi, ys, ts, x_nom, u_nom, u_hat, b, options, nvp)
%
% Multipliers adjustment method.

arguments
    f
    L
    F
    psi
    ys
    ts
    x_nom
    u_nom
    u_hat
    b
    options
    nvp.beta1 = [];
    nvp.beta2 = [];
    nvp.g = [];
    nvp.gb = [];
    nvp.hessMinSolver = [];
    nvp.a = [];
    nvp.Vx = [];
    nvp.Vxx = [];
end


nx = length(x_nom(1,:));
npsi = length(b);

% Retrieve Vb(t_0) and Vbb(t_0)
Vb0 = ys(1, 1:length(b));
Vb0 = Vb0(:);
Vbb0 = ys(1, end-length(b)^2+1:end);
Vbb0 = reshape(Vbb0, [length(b) length(b)]);
VbbVb0 = Vbb0 \ Vb0;
% Retrieve Vxb
if options.compTrick
    Vxb = ys(:, (npsi+1):(npsi+npsi*nx));
end
% Nominal endpoint error
psi_nom = full( psi.fun(x_nom(end,:), ts(end)) );

% Initialize flags
test2_enable = true;
success = false;

% Loop until satisfactory db is found
while true
    % (Re-)initialize epsilon
    epsilon = 2;

    for n = 1:10
        % Reduce db
        epsilon = epsilon/2;
        db = - epsilon * VbbVb0;

        % New trial trajectories
        if options.compTrick
            % Obtain the new control trajectory with the computational
            % trick
            [~, x_nom_new, u_nom_new] = ...
                sys_fwd(f, x_nom(1,:), ts, x_nom, u_hat, options, ...
                "db", db, "g", nvp.g, "gb", nvp.gb, ...
                "hessMinSolver", nvp.hessMinSolver, ...
                "Vx", nvp.Vx, "Vxx", nvp.Vxx, "Vxb", Vxb);
        else
            % Obtain the new control trajectory with beta1 and beta2
            [~, x_nom_new, u_nom_new] = ...
                sys_fwd(f, x_nom(1,:), ts, x_nom, u_hat, options, ...
                "db", db, "g", nvp.g, "gb", nvp.gb, ...
                "beta1", nvp.beta1, "beta2", nvp.beta2, "db", db);
        end

        % Test 1/2
        if full( norm( psi.fun(x_nom_new(end,:), ts(end) ) ) ) > norm( psi_nom )
            % Test 1 failed; reject db
            continue
        end

        if test2_enable
            % Test 2/2

            % Actual improvement in cost
            V_nom = totalCost(x_nom, u_nom, ts, L, F, b, psi);
            V = totalCost(x_nom_new, u_nom_new, ts, L, F, b + db, psi);
            DV = V - V_nom;

            % Predicted improvement in cost
            DVtilde = nvp.a(1) - (epsilon - 0.5 * epsilon^2) * Vb0.' * VbbVb0;

            gamma1 = 0.8;
            gamma2 = 1.2;
            if DV/DVtilde > gamma1 && DV/DVtilde < gamma2
                % Both tests passed; Accept db
                success = true;
                break
            else
                % Test 2 failed; reject db
                continue
            end
        else
            % Test 1 passed and test 2 disabled; accept db
            success = true;
            break
        end
    end

    if test2_enable && ~success
        % Try again with test 2 disabled
        test2_enable = false;
        continue
    elseif success
        x_nom = x_nom_new;
        u_nom = u_nom_new;
        b = b + db;
        break
    else
        warning("Multipliers adjustment method failed.")
        break
    end
end

end