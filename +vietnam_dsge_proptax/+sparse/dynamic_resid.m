function [residual, T_order, T] = dynamic_resid(y, x, params, steady_state, T_order, T)
if nargin < 6
    T_order = -1;
    T = NaN(8, 1);
end
[T_order, T] = vietnam_dsge_proptax.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
residual = NaN(18, 1);
    residual(1) = (1/y(20)) - (params(1)*(1+(1-params(23))*params(4)*y(37)/y(22)-params(5))/y(38));
    residual(2) = (1/y(20)) - (params(1)*(1+y(34))/y(38));
    residual(3) = ((1-params(23))*(1-params(4))*y(31)) - (y(20)*params(3)*T(1));
    residual(4) = (y(19)) - (T(2)^(params(6)/(params(6)-1)));
    residual(5) = (y(31)) - (T(3)*T(4));
    residual(6) = (y(28)) - (1-T(5));
    residual(7) = (y(29)) - (T(7)^(params(11)/(params(11)-1)));
    residual(8) = (y(26)/(steady_state(8))) - (T(8)/y(32));
    residual(9) = (y(23)) - (y(26)+y(5)*(1-params(12)));
    residual(10) = (y(30)) - (y(12)*(1+params(19)*params(20)*(params(14)-y(28))/params(14)));
    residual(11) = (y(27)/(steady_state(9))) - (T(8)*exp(y(36)));
    residual(12) = (y(24)) - (y(27)+y(6)*(1-params(13)));
    residual(13) = (y(22)) - (y(4)*(1-params(5))+y(25));
    residual(14) = (y(33)) - (y(19)+(1+y(16))*y(15)-y(20)-y(25)-y(26)*y(32)-y(27));
    residual(15) = (y(34)) - (params(25)+params(26)*(exp(params(27)-y(33))-1));
    residual(16) = (log(y(32))) - (params(21)*log(y(14))+x(2));
    residual(17) = (y(35)) - (params(17)*y(17)+x(1));
    residual(18) = (y(36)) - (params(18)*y(18)+x(3));
end
