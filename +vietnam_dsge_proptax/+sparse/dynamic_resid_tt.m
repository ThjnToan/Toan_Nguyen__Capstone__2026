function [T_order, T] = dynamic_resid_tt(y, x, params, steady_state, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 8
    T = [T; NaN(8 - size(T, 1), 1)];
end
T(1) = y(21)^(1+params(2));
T(2) = (1-params(7))*y(31)^((params(6)-1)/params(6))+params(7)*params(28)^((params(6)-1)/params(6));
T(3) = params(30)*(y(28)*y(4))^params(4);
T(4) = y(21)^(1-params(4));
T(5) = exp((-params(8))*y(29)/(params(9)*params(29)*exp(y(35))));
T(6) = (params(11)-1)/params(11);
T(7) = params(10)*(y(12)*y(5))^T(6)+(1-params(10))*y(6)^T(6);
T(8) = (params(14)/y(28))^params(24);
end
