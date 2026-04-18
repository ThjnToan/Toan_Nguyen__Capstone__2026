function g1 = static_g1(T, y, x, params, T_flag)
% function g1 = static_g1(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g1
%

if T_flag
    T = vietnam_dsge_proptax.static_g1_tt(T, y, x, params);
end
g1 = zeros(18, 18);
g1(1,1)=(-(params(1)*(1-params(23))*params(4)/y(4)/y(2)));
g1(1,2)=(-1)/(y(2)*y(2))-(-(params(1)*(1+(1-params(23))*params(4)*y(1)/y(4)-params(5))))/(y(2)*y(2));
g1(1,4)=(-(params(1)*(-((1-params(23))*params(4)*y(1)))/(y(4)*y(4))/y(2)));
g1(2,2)=(-1)/(y(2)*y(2))-(-(params(1)*(1+y(16))))/(y(2)*y(2));
g1(2,16)=(-(params(1)/y(2)));
g1(3,2)=(-(params(3)*T(1)));
g1(3,3)=(-(y(2)*params(3)*getPowerDeriv(y(3),1+params(2),1)));
g1(3,13)=(1-params(23))*(1-params(4));
g1(4,1)=1;
g1(4,13)=(-((1-params(7))*getPowerDeriv(y(13),(params(6)-1)/params(6),1)*getPowerDeriv(T(2),params(6)/(params(6)-1),1)));
g1(5,3)=(-(T(3)*getPowerDeriv(y(3),1-params(4),1)));
g1(5,4)=(-(T(4)*params(30)*y(10)*T(9)));
g1(5,10)=(-(T(4)*params(30)*y(4)*T(9)));
g1(5,13)=1;
g1(6,10)=1;
g1(6,11)=T(5)*(-params(8))/(params(9)*params(29)*exp(y(17)));
g1(6,17)=T(5)*(-((-params(8))*y(11)*params(9)*params(29)*exp(y(17))))/(params(9)*params(29)*exp(y(17))*params(9)*params(29)*exp(y(17)));
g1(7,5)=(-(params(10)*y(12)*T(10)*T(11)));
g1(7,6)=(-(T(11)*(1-params(10))*getPowerDeriv(y(6),T(6),1)));
g1(7,11)=1;
g1(7,12)=(-(T(11)*params(10)*y(5)*T(10)));
g1(8,8)=((y(8))-y(8))/((y(8))*(y(8)));
g1(8,10)=(-(T(12)/y(14)));
g1(8,14)=(-((-T(8))/(y(14)*y(14))));
g1(9,5)=1-(1-params(12));
g1(9,8)=(-1);
g1(10,10)=(-(y(12)*(-(params(19)*params(20)))/params(14)));
g1(10,12)=1-(1+params(19)*params(20)*(params(14)-y(10))/params(14));
g1(11,9)=((y(9))-y(9))/((y(9))*(y(9)));
g1(11,10)=(-(exp(y(18))*T(12)));
g1(11,18)=(-(T(8)*exp(y(18))));
g1(12,6)=1-(1-params(13));
g1(12,9)=(-1);
g1(13,4)=1-(1-params(5));
g1(13,7)=(-1);
g1(14,1)=(-1);
g1(14,2)=1;
g1(14,7)=1;
g1(14,8)=y(14);
g1(14,9)=1;
g1(14,14)=y(8);
g1(14,15)=1-(1+y(16));
g1(14,16)=(-y(15));
g1(15,15)=(-(params(26)*(-exp(params(27)-y(15)))));
g1(15,16)=1;
g1(16,14)=1/y(14)-params(21)*1/y(14);
g1(17,17)=1-params(17);
g1(18,18)=1-params(18);

end
