function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
% function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = vietnam_dsge_proptax.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(18, 1);
    residual(1) = (1/y(11)) - (params(1)*(1+(1-params(23))*params(4)*y(28)/y(13)-params(5))/y(29));
    residual(2) = (1/y(11)) - (params(1)*(1+y(25))/y(29));
    residual(3) = ((1-params(23))*(1-params(4))*y(22)) - (y(11)*params(3)*T(1));
    residual(4) = (y(10)) - (T(2)^(params(6)/(params(6)-1)));
    residual(5) = (y(22)) - (T(3)*T(4));
    residual(6) = (y(19)) - (1-T(5));
    residual(7) = (y(20)) - (T(7)^(params(11)/(params(11)-1)));
    residual(8) = (y(17)/(steady_state(8))) - (T(8)/y(23));
    residual(9) = (y(14)) - (y(17)+y(2)*(1-params(12)));
    residual(10) = (y(21)) - (y(4)*(1+params(19)*params(20)*(params(14)-y(19))/params(14)));
    residual(11) = (y(18)/(steady_state(9))) - (T(8)*exp(y(27)));
    residual(12) = (y(15)) - (y(18)+y(3)*(1-params(13)));
    residual(13) = (y(13)) - (y(1)*(1-params(5))+y(16));
    residual(14) = (y(24)) - (y(10)+(1+y(7))*y(6)-y(11)-y(16)-y(17)*y(23)-y(18));
    residual(15) = (y(25)) - (params(25)+params(26)*(exp(params(27)-y(24))-1));
    residual(16) = (log(y(23))) - (params(21)*log(y(5))+x(it_, 2));
    residual(17) = (y(26)) - (params(17)*y(8)+x(it_, 1));
    residual(18) = (y(27)) - (params(18)*y(9)+x(it_, 3));

end
