function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
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
%   g1
%

if T_flag
    T = vietnam_dsge_proptax.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(18, 32);
g1(1,28)=(-(params(1)*(1-params(23))*params(4)/y(13)/y(29)));
g1(1,11)=(-1)/(y(11)*y(11));
g1(1,29)=(-((-(params(1)*(1+(1-params(23))*params(4)*y(28)/y(13)-params(5))))/(y(29)*y(29))));
g1(1,13)=(-(params(1)*(-((1-params(23))*params(4)*y(28)))/(y(13)*y(13))/y(29)));
g1(2,11)=(-1)/(y(11)*y(11));
g1(2,29)=(-((-(params(1)*(1+y(25))))/(y(29)*y(29))));
g1(2,25)=(-(params(1)/y(29)));
g1(3,11)=(-(params(3)*T(1)));
g1(3,12)=(-(y(11)*params(3)*getPowerDeriv(y(12),1+params(2),1)));
g1(3,22)=(1-params(23))*(1-params(4));
g1(4,10)=1;
g1(4,22)=(-((1-params(7))*getPowerDeriv(y(22),(params(6)-1)/params(6),1)*getPowerDeriv(T(2),params(6)/(params(6)-1),1)));
g1(5,12)=(-(T(3)*getPowerDeriv(y(12),1-params(4),1)));
g1(5,1)=(-(T(4)*params(30)*y(19)*T(9)));
g1(5,19)=(-(T(4)*params(30)*y(1)*T(9)));
g1(5,22)=1;
g1(6,19)=1;
g1(6,20)=T(5)*(-params(8))/(params(9)*params(29)*exp(y(26)));
g1(6,26)=T(5)*(-((-params(8))*y(20)*params(9)*params(29)*exp(y(26))))/(params(9)*params(29)*exp(y(26))*params(9)*params(29)*exp(y(26)));
g1(7,2)=(-(params(10)*y(4)*T(10)*T(11)));
g1(7,3)=(-(T(11)*(1-params(10))*getPowerDeriv(y(3),T(6),1)));
g1(7,20)=1;
g1(7,4)=(-(T(11)*params(10)*y(2)*T(10)));
g1(8,17)=1/(steady_state(8));
g1(8,19)=(-(T(12)/y(23)));
g1(8,23)=(-((-T(8))/(y(23)*y(23))));
g1(9,2)=(-(1-params(12)));
g1(9,14)=1;
g1(9,17)=(-1);
g1(10,19)=(-(y(4)*(-(params(19)*params(20)))/params(14)));
g1(10,4)=(-(1+params(19)*params(20)*(params(14)-y(19))/params(14)));
g1(10,21)=1;
g1(11,18)=1/(steady_state(9));
g1(11,19)=(-(exp(y(27))*T(12)));
g1(11,27)=(-(T(8)*exp(y(27))));
g1(12,3)=(-(1-params(13)));
g1(12,15)=1;
g1(12,18)=(-1);
g1(13,1)=(-(1-params(5)));
g1(13,13)=1;
g1(13,16)=(-1);
g1(14,10)=(-1);
g1(14,11)=1;
g1(14,16)=1;
g1(14,17)=y(23);
g1(14,18)=1;
g1(14,23)=y(17);
g1(14,6)=(-(1+y(7)));
g1(14,24)=1;
g1(14,7)=(-y(6));
g1(15,24)=(-(params(26)*(-exp(params(27)-y(24)))));
g1(15,25)=1;
g1(16,5)=(-(params(21)*1/y(5)));
g1(16,23)=1/y(23);
g1(16,31)=(-1);
g1(17,8)=(-params(17));
g1(17,26)=1;
g1(17,30)=(-1);
g1(18,9)=(-params(18));
g1(18,27)=1;
g1(18,32)=(-1);

end
