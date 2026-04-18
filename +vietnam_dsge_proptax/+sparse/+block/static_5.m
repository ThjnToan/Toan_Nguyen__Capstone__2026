function [y, T, residual, g1] = static_5(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(8, 1);
  T(10)=(1-params(7))*y(13)^((params(6)-1)/params(6))+params(7)*params(28)^((params(6)-1)/params(6));
  residual(1)=(y(1))-(T(10)^(params(6)/(params(6)-1)));
  T(11)=params(30)*(y(4)*y(10))^params(4);
  T(12)=y(3)^(1-params(4));
  residual(2)=(y(13))-(T(11)*T(12));
  residual(3)=(y(4))-(y(4)*(1-params(5))+y(7));
  residual(4)=(y(15))-(y(1)+(1+y(16))*y(15)-y(2)-y(7)-y(8)*y(14)-y(9));
  T(13)=exp(params(27)-y(15));
  residual(5)=(y(16))-(params(25)+params(26)*(T(13)-1));
  residual(6)=(1/y(2))-(params(1)*(1+(1-params(23))*params(4)*y(1)/y(4)-params(5))/y(2));
  residual(7)=(1/y(2))-(params(1)*(1+y(16))/y(2));
  T(14)=y(3)^(1+params(2));
  residual(8)=((1-params(23))*(1-params(4))*y(13))-(y(2)*params(3)*T(14));
if nargout > 3
    g1_v = NaN(22, 1);
g1_v(1)=1;
g1_v(2)=(-1);
g1_v(3)=(-(params(1)*(1-params(23))*params(4)/y(4)/y(2)));
g1_v(4)=(-((1-params(7))*getPowerDeriv(y(13),(params(6)-1)/params(6),1)*getPowerDeriv(T(10),params(6)/(params(6)-1),1)));
g1_v(5)=1;
g1_v(6)=(1-params(23))*(1-params(4));
g1_v(7)=(-1);
g1_v(8)=1;
g1_v(9)=1;
g1_v(10)=(-1)/(y(2)*y(2))-(-(params(1)*(1+(1-params(23))*params(4)*y(1)/y(4)-params(5))))/(y(2)*y(2));
g1_v(11)=(-1)/(y(2)*y(2))-(-(params(1)*(1+y(16))))/(y(2)*y(2));
g1_v(12)=(-(params(3)*T(14)));
g1_v(13)=1-(1+y(16));
g1_v(14)=(-(params(26)*(-T(13))));
g1_v(15)=(-(T(12)*params(30)*y(10)*getPowerDeriv(y(4)*y(10),params(4),1)));
g1_v(16)=1-(1-params(5));
g1_v(17)=(-(params(1)*(-((1-params(23))*params(4)*y(1)))/(y(4)*y(4))/y(2)));
g1_v(18)=(-y(15));
g1_v(19)=1;
g1_v(20)=(-(params(1)/y(2)));
g1_v(21)=(-(T(11)*getPowerDeriv(y(3),1-params(4),1)));
g1_v(22)=(-(y(2)*params(3)*getPowerDeriv(y(3),1+params(2),1)));
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 8, 8);
end
end
