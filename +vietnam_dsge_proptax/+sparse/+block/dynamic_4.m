function [y, T, residual, g1] = dynamic_4(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(8, 1);
  T(9)=params(30)*(y(28)*y(4))^params(4);
  T(10)=y(21)^(1-params(4));
  residual(1)=(y(31))-(T(9)*T(10));
  residual(2)=(y(33))-(y(19)+(1+y(16))*y(15)-y(20)-y(25)-y(26)*y(32)-y(27));
  T(11)=y(21)^(1+params(2));
  residual(3)=((1-params(23))*(1-params(4))*y(31))-(y(20)*params(3)*T(11));
  residual(4)=(y(22))-(y(4)*(1-params(5))+y(25));
  T(12)=exp(params(27)-y(33));
  residual(5)=(y(34))-(params(25)+params(26)*(T(12)-1));
  residual(6)=(1/y(20))-(params(1)*(1+y(34))/y(38));
  T(13)=(1-params(7))*y(31)^((params(6)-1)/params(6))+params(7)*params(28)^((params(6)-1)/params(6));
  residual(7)=(y(19))-(T(13)^(params(6)/(params(6)-1)));
  residual(8)=(1/y(20))-(params(1)*(1+(1-params(23))*params(4)*y(37)/y(22)-params(5))/y(38));
if nargout > 3
    g1_v = NaN(26, 1);
g1_v(1)=(-(T(10)*params(30)*y(28)*getPowerDeriv(y(28)*y(4),params(4),1)));
g1_v(2)=(-(1-params(5)));
g1_v(3)=(-(1+y(16)));
g1_v(4)=(-y(15));
g1_v(5)=1;
g1_v(6)=(1-params(23))*(1-params(4));
g1_v(7)=(-((1-params(7))*getPowerDeriv(y(31),(params(6)-1)/params(6),1)*getPowerDeriv(T(13),params(6)/(params(6)-1),1)));
g1_v(8)=1;
g1_v(9)=(-1);
g1_v(10)=(-(T(9)*getPowerDeriv(y(21),1-params(4),1)));
g1_v(11)=(-(y(20)*params(3)*getPowerDeriv(y(21),1+params(2),1)));
g1_v(12)=1;
g1_v(13)=(-(params(1)*(-((1-params(23))*params(4)*y(37)))/(y(22)*y(22))/y(38)));
g1_v(14)=1;
g1_v(15)=(-(params(26)*(-T(12))));
g1_v(16)=1;
g1_v(17)=(-(params(1)/y(38)));
g1_v(18)=(-1);
g1_v(19)=1;
g1_v(20)=1;
g1_v(21)=(-(params(3)*T(11)));
g1_v(22)=(-1)/(y(20)*y(20));
g1_v(23)=(-1)/(y(20)*y(20));
g1_v(24)=(-(params(1)*(1-params(23))*params(4)/y(22)/y(38)));
g1_v(25)=(-((-(params(1)*(1+y(34))))/(y(38)*y(38))));
g1_v(26)=(-((-(params(1)*(1+(1-params(23))*params(4)*y(37)/y(22)-params(5))))/(y(38)*y(38))));
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 8, 24);
end
end
