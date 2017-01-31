function  res = TV_Spat()


% spatial finite-differencing operator along x and y dimension
%

res.adjoint = 0;
res = class(res,'TV_Spat');

