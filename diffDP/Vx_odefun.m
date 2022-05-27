function [dydt, beta1, Z, A] = Vx_odefun(t, y, x_nom, u_nom, u_hat, f, H, g_hat, nx, nu, options)
%Vx_odefun(t, y, x_nom, u_nom, u_hat, f, H, g_hat, nx, nu, options)
% 
% Base equations for a, Vx, Vxx.

% Check whether the solution is unbounded
if any(isnan(y))
    dydt = nan(size(y));
    beta1 = nan(nu,nx);
    Z = nan(nu);
    A = nan(nu);
    return
end

x_nom = x_nom(:);
u_nom = u_nom(:);
Vx = reshape(y(2:(1+nx)), [1 nx]);
Vxx = reshape(y(2+nx:end), [nx nx]);

% Ingredients
f_hat = full( f.fun(x_nom, u_hat, t) );
H_hat = full( H.fun(x_nom, u_hat, Vx, t) );

f_bar = full( f.fun(x_nom, u_nom, t) );
H_bar = full( H.fun(x_nom, u_nom, Vx, t) );

fx = full( f.fx(x_nom, u_hat, t) );
fu = full( f.fu(x_nom, u_hat, t) );

Hx = full( H.Hx(x_nom, u_hat, Vx, t) )';
if options.nullHxx
    Hxx = zeros(nx);
else
    Hxx = full( H.Hxx(x_nom, u_hat, Vx, t) );
end
Hux = full( H.Hux(x_nom, u_hat, Vx, t) );
Huu = full( H.Huu(x_nom, u_hat, Vx, t) );

p_hat = length(g_hat); % number of active constraints

if p_hat > 0
    gx_hat = zeros(p_hat, nx);
    gu_hat = zeros(p_hat, nu);
    for n = 1:p_hat
        gx_hat(n,:) = full( g_hat{n}.gx(x_nom, u_hat, t) );
        gu_hat(n,:) = full( g_hat{n}.gu(x_nom, u_hat, t) );
    end

    Hu = full( H.Hu(x_nom, u_nom, Vx, t) )';

    % Multiplier function
    mu = - inv(gu_hat * gu_hat.') * gu_hat * Hu;

    % Evaluate guu_hat
    gxx_hat = zeros(p_hat, nx, nx);
    guu_hat = zeros(p_hat, nu, nu);
    gux_hat = zeros(p_hat, nu, nx);
    for n = 1:p_hat
        gxx_hat(n,:,:) = full( g_hat{n}.gxx(x_nom, u_hat, t) );
        gux_hat(n,:,:) = full( g_hat{n}.gux(x_nom, u_hat, t) );
        guu_hat(n,:,:) = full( g_hat{n}.guu(x_nom, u_hat, t) );
    end

    % Evaluate Q and Z
    A = inv(Huu +  mu*guu_hat);
    Q = inv(gu_hat * A *  gu_hat.') * gu_hat * A;
    Z = eye(nu) - gu_hat.' * Q ;

    mu_gxx = tensor_times_vector(gxx_hat, mu, 1);
    mu_guu = tensor_times_vector(guu_hat, mu, 1);
    mu_gux = tensor_times_vector(gux_hat, mu, 1);

    beta1 = - A * Z * (Hux + mu_gux + fu.'*Vxx) - Q .' * gx_hat;
else
    A = inv(Huu);
    Q = 0;
    Z = eye(nu);

    mu = 0;
    gx_hat = 0;
    mu_gxx = 0;
    mu_guu = 0;
    mu_gux = 0;

    beta1 = - (Z * A) * (Hux + fu.'*Vxx) - Q*gx_hat;

end

% Differential equations da/dt, dVx/dt, dVxx/dt
a_dot = H_hat - H_bar;
Vx_dot = Hx + Vxx * ( f_hat - f_bar ) + (gx_hat.' * mu);
Vxx_dot = Hxx + mu_gxx + fx.' * Vxx + Vxx * fx + ...
    + beta1.' * (Huu + mu_guu) * beta1 + ...
    + beta1.' * (Hux + mu_gux + fu.'*Vxx) + ...
    + (Hux + mu_gux + fu.'*Vxx).' * beta1;

dydt = [a_dot; Vx_dot(:); Vxx_dot(:)];
dydt = - dydt;

end
