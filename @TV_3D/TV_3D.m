function  res = TV_3D(stepsize)
if nargin < 1
    stepsize = [1,1,1];
end

% spatial and temporal finite-differencing operator along x, y and t dimension
%
% [res.D,res.Dt]=defDDt(stepsize);
res.stepsize = stepsize;
res.adjoint = 0;
res = class(res,'TV_3D');
