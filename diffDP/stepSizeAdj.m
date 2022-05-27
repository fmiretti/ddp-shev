function [x_nom, u_nom, halt] = stepSizeAdj(f, L, F, psi, ys, ts, x_nom, u_nom, u_hat, b, N_eff, options, nvp)
%stepSizeAdj
% 
% Step-size adjustment method.
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
    N_eff
    options
    nvp.beta1 = [];
    nvp.g = [];
    nvp.gb = [];
    nvp.hessMinSolver = [];
end

nx = size(x_nom,2);

C = 0.5;
to_next_iter = false;
halt = false;
while true
    N1 = int32(2 - N_eff);
    while N1 < (N_eff-1)
        % Adjust N1
        N1 = ceil( (N_eff - N1) / 2 + N1 );
        % Skip unbouded segments
        if any(isnan(u_hat( N1, : )))
            continue
        end
        
        % (Re-)evaluate new nominal trajectories
        ts_t0_t1 = ts( 1:N1-1 )'; % Time nodes from t0 to t1
        ts_t1_tf = ts( N1:end )'; % Time nodes from t1 to tf

        x_nom_t0_t1 = x_nom( 1:N1-1, : );
        u_nom_t0_t1 = u_nom( 1:N1-1, : );
        x_nom_t1_tf = x_nom( N1:end, : );
        u_hat_t1_tf = u_hat( N1:end, : );
        
        % Apply the new control trajectory from t1 to the end
        if options.compTrick
            % Obtain the new control trajectory with the computational
            % trick
            Vx = ys( N1:end, 2:nx+1 );
            Vxx = ys( N1:end, nx+2:end );
            [ts_t1_tf, x_nom_new_t1_tf, u_nom_new_t1_tf] = ...
                sys_fwd(f, x_nom(N1,:), ts_t1_tf, x_nom_t1_tf, u_hat_t1_tf, options, ...
                "g", nvp.g, "gb", nvp.gb, ...
                "hessMinSolver", nvp.hessMinSolver, "Vx", Vx, "Vxx", Vxx);
        else
            % Obtain the new control trajectory with beta1 * dx
            beta1_t1_tf = nvp.beta1( N1:end, : );
            [ts_t1_tf, x_nom_new_t1_tf, u_nom_new_t1_tf] = ...
                sys_fwd(f, x_nom(N1,:), ts_t1_tf, x_nom_t1_tf, u_hat_t1_tf, options, ...
                "beta1", beta1_t1_tf, "g", nvp.g, "gb", nvp.gb);
        end

        x_nom_new = [x_nom_t0_t1; x_nom_new_t1_tf];
        u_nom_new = [u_nom_t0_t1; u_nom_new_t1_tf];

        % Actual improvement in cost
        V_nom = totalCost(x_nom_t1_tf, u_nom( N1:end, : ), ts_t1_tf, L, F, b, psi);
        V = totalCost(x_nom_new_t1_tf, u_nom_new_t1_tf, ts_t1_tf, L, F, b, psi);
        DV1 = V_nom - V;

        % Predicted improvement in cost
        a1 = abs( ys(N1,1) );

        % If the solution became unbounded, V0_new will be nan. The
        % step size adjustment method must be used.
        if isnan(DV1)
            DV1 = 0;
        end
        if isnan(a1)
            a1 = 1;
        end

        % Check improvement
        if DV1/a1 > C
            % Good improvement in cost achieved. Move to the next iter.
            to_next_iter = true;
            break
        end
    end

    %% Flow control
    if to_next_iter
        % Good improvement in cost achieved. Move to the next iter.
        % Save new trajectories and cost
        x_nom = x_nom_new;
        u_nom = u_nom_new;
        return
    end

    % No satisfactory t1 found.
    if C == 0
        % No improvement in trajectory attainable.
        warning("Unable to solve the free endpoint problem.")
        halt = true;
        return
    else
        % Retry the step-size adjustment method with less stringent
        % acceptance criterion
        C = 0;
    end
end

end