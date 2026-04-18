function [residual, T_order, T] = static_resid(y, x, params, T_order, T)
if nargin < 5
    T_order = -1;
    T = NaN(8, 1);
end
[T_order, T] = vietnam_dsge_proptax.sparse.static_resid_tt(y, x, params, T_order, T);
residual = NaN(18, 1);
    residual(1) = (1/y(2)) - (params(1)*(1+(1-params(23))*params(4)*y(1)/y(4)-params(5))/y(2));
    residual(2) = (1/y(2)) - (params(1)*(1+y(16))/y(2));
    residual(3) = ((1-params(23))*(1-params(4))*y(13)) - (y(2)*params(3)*T(1));
    residual(4) = (y(1)) - (T(2)^(params(6)/(params(6)-1)));
    residual(5) = (y(13)) - (T(3)*T(4));
    residual(6) = (y(10)) - (1-T(5));
    residual(7) = (y(11)) - (T(7)^(params(11)/(params(11)-1)));
    residual(8) = (y(8)/(y(8))) - (T(8)/y(14));
    residual(9) = (y(5)) - (y(8)+y(5)*(1-params(12)));
    residual(10) = (y(12)) - (y(12)*(1+params(19)*params(20)*(params(14)-y(10))/params(14)));
    residual(11) = (y(9)/(y(9))) - (T(8)*exp(y(18)));
    residual(12) = (y(6)) - (y(9)+y(6)*(1-params(13)));
    residual(13) = (y(4)) - (y(4)*(1-params(5))+y(7));
    residual(14) = (y(15)) - (y(1)+(1+y(16))*y(15)-y(2)-y(7)-y(8)*y(14)-y(9));
    residual(15) = (y(16)) - (params(25)+params(26)*(exp(params(27)-y(15))-1));
    residual(16) = (log(y(14))) - (log(y(14))*params(21)+x(2));
    residual(17) = (y(17)) - (y(17)*params(17)+x(1));
    residual(18) = (y(18)) - (y(18)*params(18)+x(3));
end
