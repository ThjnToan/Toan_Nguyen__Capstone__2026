function [y, T, residual, g1] = dynamic_3(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(6, 1);
  T(1)=params(9)*params(29)*exp(y(35));
  T(2)=exp((-params(8))*y(29)/T(1));
  y(28)=1-T(2);
  T(3)=params(10)*(y(12)*y(5))^((params(11)-1)/params(11))+(1-params(10))*y(6)^((params(11)-1)/params(11));
  residual(1)=(y(29))-(T(3)^(params(11)/(params(11)-1)));
  T(4)=(params(14)/y(28))^params(24);
  residual(2)=(y(26)/(steady_state(8)))-(T(4)/y(32));
  T(5)=exp(y(36));
  residual(3)=(y(27)/(steady_state(9)))-(T(4)*T(5));
  residual(4)=(y(23))-(y(26)+y(5)*(1-params(12)));
  residual(5)=(y(30))-(y(12)*(1+params(19)*params(20)*(params(14)-y(28))/params(14)));
  residual(6)=(y(24))-(y(27)+y(6)*(1-params(13)));
  T(6)=getPowerDeriv(y(12)*y(5),(params(11)-1)/params(11),1);
  T(7)=getPowerDeriv(T(3),params(11)/(params(11)-1),1);
  T(8)=getPowerDeriv(params(14)/y(28),params(24),1);
if nargout > 3
    g1_v = NaN(11, 1);
g1_v(1)=1;
g1_v(2)=(-(T(8)*(-(params(14)*(-(T(2)*(-params(8))/T(1)))))/(y(28)*y(28))/y(32)));
g1_v(3)=(-(T(5)*T(8)*(-(params(14)*(-(T(2)*(-params(8))/T(1)))))/(y(28)*y(28))));
g1_v(4)=(-(y(12)*params(19)*params(20)*T(2)*(-params(8))/T(1)/params(14)));
g1_v(5)=1/(steady_state(8));
g1_v(6)=(-1);
g1_v(7)=1/(steady_state(9));
g1_v(8)=(-1);
g1_v(9)=1;
g1_v(10)=1;
g1_v(11)=1;
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 6, 6);
end
end
