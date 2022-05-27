function [ts, y, beta2] = Vb_bwd(ts, y0, x_nom, u_nom, u_hat, beta1, A, Z, f, nx, nu, npsi, options )
%Vb_bwd
%
% Integrate base equations for Vb, Vxb, Vbb.

% storage space for beta2
beta2 = zeros(length(ts), nu, npsi);

% Storage space for y = [Vb(:), Vxb(:), Vbb(:)];
y = zeros(length(ts), length(y0));

% Initial conditions
y(end, :) = y0;

for n = length(ts):-1:1
    yy = y(n, :);
    yy = yy(:);

    xx_nom =  x_nom(n, :);
    uu_nom = u_nom(n, :);
    uu_hat = u_hat(n, :);
    bbeta1 = beta1(n, :);
    AA = A(n, :);
    ZZ = Z(n, :);

    [k_1, beta2(n,:,:)] = Vb_odefun( ts(n), yy, xx_nom, uu_nom, uu_hat, bbeta1, AA, ZZ, f, nx, nu, npsi );

    if n>1
        % Timestep
        h = ts(n-1) - ts(n);
        switch options.integrator
            case "rk4"
                k_2 = Vb_odefun( ts(n) + h/2, yy + k_1*h/2, xx_nom, uu_nom, uu_hat, bbeta1, AA, ZZ, f, nx, nu, npsi );
                k_3 = Vb_odefun( ts(n) + h/2, yy + k_2*h/2, xx_nom, uu_nom, uu_hat, bbeta1, AA, ZZ, f, nx, nu, npsi );
                k_4 = Vb_odefun( ts(n) + h, yy + k_3*h, xx_nom, uu_nom, uu_hat, bbeta1, AA, ZZ, f, nx, nu, npsi );
                y(n-1, :) = yy + (1/6)*( k_1 + 2*k_2 + 2*k_3 + k_4 )*h;
            case "fwdEuler"
                y(n-1, :) = yy + k_1*h;
        end
    end
end

end