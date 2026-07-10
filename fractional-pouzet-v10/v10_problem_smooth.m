function out = v10_problem_smooth(t,y,flag)
% Smooth manufactured scalar Caputo problem used in the homogeneous benchmark.
% Exact solution: y(t)=1+(t/T)^(5+alpha).
% This file is vectorized as required by fhbvm.
global V10_ALPHA V10_T
if nargin==0
    out = V10_ALPHA;
    return
end
alpha = V10_ALPHA; T = V10_T;
if nargin<3
    tt = t(:);
    ye = 1 + (tt./T).^(5+alpha);
    cap = gamma(6+alpha)/gamma(6) .* tt.^5 ./ T.^(5+alpha);
    q = cap + ye.^2;
    if isscalar(y)
        out = -y.^2 + q(1);
    else
        out = -y.^2 + q;
    end
else
    out = -2*y;
end
end
