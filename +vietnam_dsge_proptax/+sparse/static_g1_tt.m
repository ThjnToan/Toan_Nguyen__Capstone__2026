function [T_order, T] = static_g1_tt(y, x, params, T_order, T)
if T_order >= 1
    return
end
[T_order, T] = vietnam_dsge_proptax.sparse.static_resid_tt(y, x, params, T_order, T);
T_order = 1;
if size(T, 1) < 12
    T = [T; NaN(12 - size(T, 1), 1)];
end
T(9) = getPowerDeriv(y(4)*y(10),params(4),1);
T(10) = getPowerDeriv(y(12)*y(5),T(6),1);
T(11) = getPowerDeriv(T(7),params(11)/(params(11)-1),1);
T(12) = (-params(14))/(y(10)*y(10))*getPowerDeriv(params(14)/y(10),params(24),1);
end
