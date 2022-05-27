function [x_next, stage_cost, unfeas] = dp_model(x, u, w, f, L)
% x{1}: SOE
% u{1}: normalized genset power
% w{1}: time

dt = w{2};

% State dynamics
if length(x) == 1
    dxdt = f(x{1}, u{1}, w{1});
    x_next{1} = x{1} + dxdt .* dt;
    unfeas = false;
else
    dxdt{1} = f{1}(x, u{1}, w{1});
    dxdt{2} = f{2}(x, u{1}, w{1});
    
    unfeas = false(size(dxdt{1}));
    unfeas( imag(dxdt{1}) ~= 0 | imag(dxdt{2}) ~= 0 ) = true;
    dxdt{1} = real(dxdt{1});
    dxdt{2} = real(dxdt{2});

    x_next{1} = x{1} + dxdt{1} .* dt;
    x_next{2} = x{2} + dxdt{2} .* dt;

end

% Cost
stage_cost = L(x{1}, u{1}, w{1}) .* dt;

end