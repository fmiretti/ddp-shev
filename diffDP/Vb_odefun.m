function [dydt, beta2] = Vb_odefun(t, y, x_nom, u_nom, u_hat, beta1, A, Z, f, nx, nu, npsi)
%Vb_odefun
%
% Base equations for Vb, Vxb, Vbb.

x_nom = x_nom(:);
u_nom = u_nom(:);

Vxb = reshape(y(npsi+1:(npsi+nx*npsi)), nx, npsi);

% Ingredients for the differential equations
f_hat = full( f.fun(x_nom, u_hat, t) );

f_bar = full( f.fun(x_nom, u_nom, t) );
% Derivatives
fx = full( f.fx(x_nom, u_hat, t) );
fu = full( f.fu(x_nom, u_hat, t) );

% Beta2
beta2 = - A * Z * fu.' * Vxb;
% Differential equations da/dt, dVx/dt, dVxx/dt
Vb_dot = Vxb.' * ( f_hat - f_bar ) ;
Vxb_dot = ( fx.' + beta1.' * fu.' ) * Vxb;
if Z ~= eye(nu)
    Vbb_dot = - Vxb.' * fu * Z.' * A * Z * fu.' * Vxb;
else
    Vbb_dot = - Vxb.' * fu * A * fu.' * Vxb;
end
dydt = [Vb_dot(:); Vxb_dot(:); Vbb_dot(:)];
dydt = - dydt;

end