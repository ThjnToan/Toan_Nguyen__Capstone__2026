function [y, T, residual, g1] = static_4(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(7, 1);
  T(2)=params(9)*params(29)*exp(y(17));
  T(3)=exp((-params(8))*y(11)/T(2));
  residual(1)=(y(10))-(1-T(3));
  T(4)=params(10)*(y(12)*y(5))^((params(11)-1)/params(11))+(1-params(10))*y(6)^((params(11)-1)/params(11));
  residual(2)=(y(11))-(T(4)^(params(11)/(params(11)-1)));
  T(5)=(params(14)/y(10))^params(24);
  residual(3)=(y(8)/(y(8)))-(T(5)/y(14));
  residual(4)=(y(5))-(y(8)+y(5)*(1-params(12)));
  residual(5)=(y(12))-(y(12)*(1+params(19)*params(20)*(params(14)-y(10))/params(14)));
  T(6)=exp(y(18));
  residual(6)=(y(9)/(y(9)))-(T(5)*T(6));
  residual(7)=(y(6))-(y(9)+y(6)*(1-params(13)));
  T(7)=getPowerDeriv(y(12)*y(5),(params(11)-1)/params(11),1);
  T(8)=getPowerDeriv(T(4),params(11)/(params(11)-1),1);
  T(9)=(-params(14))/(y(10)*y(10))*getPowerDeriv(params(14)/y(10),params(24),1);
if nargout > 3
    g1_v = NaN(16, 1);
g1_v(1)=T(3)*(-params(8))/T(2);
g1_v(2)=1;
g1_v(3)=(-(T(8)*params(10)*y(5)*T(7)));
g1_v(4)=1-(1+params(19)*params(20)*(params(14)-y(10))/params(14));
g1_v(5)=((y(8))-y(8))/((y(8))*(y(8)));
g1_v(6)=(-1);
g1_v(7)=(-(params(10)*y(12)*T(7)*T(8)));
g1_v(8)=1-(1-params(12));
g1_v(9)=1;
g1_v(10)=(-(T(9)/y(14)));
g1_v(11)=(-(y(12)*(-(params(19)*params(20)))/params(14)));
g1_v(12)=(-(T(6)*T(9)));
g1_v(13)=((y(9))-y(9))/((y(9))*(y(9)));
g1_v(14)=(-1);
g1_v(15)=(-(T(8)*(1-params(10))*getPowerDeriv(y(6),(params(11)-1)/params(11),1)));
g1_v(16)=1-(1-params(13));
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 7, 7);
end
end
