function [y, T] = dynamic_2(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(35)=params(17)*y(17)+x(1);
  y(36)=params(18)*y(18)+x(3);
end
