function V0 = totalCost(x, u, ts, L, F, b, psi)
%totalCost
%
% Evaluate total cost.
%
% Notes
% -----
% To evaluate the true cost set b to 0.

if sum(b) == 0 || isempty(b)
    b = zeros(size( full( psi.fun(x(1,:), 0) ) ));
end

Ls = zeros(length(ts), 1);
for n = 1:length(ts)
    Ls(n) = full( L( x(n,:), u(n,:), ts(n) ) );
end
V0 = trapz(ts, Ls) + full( F.fun( x(end,:), ts(end) ) ) + b.' * full( psi.fun(x(end,:), ts(end)) );

end