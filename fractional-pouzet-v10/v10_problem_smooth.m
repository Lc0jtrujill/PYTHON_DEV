function out = v10_problem_smooth(t,y,flag)
% Smooth manufactured scalar Caputo problem used in the homogeneous benchmark.
% Exact solution: y(t)=1+(t/T)^(5+alpha).
% The start y=1 has error O(h^(5+alpha)), above the P2--P4 target orders.
global V10_ALPHA V10_T
if nargin==0
    out = V10_ALPHA;
    return
end
alpha = V10_ALPHA; T = V10_T;
tt = t(:);
ye = 1 + (tt./T).^(5+alpha);
cap = gamma(6+alpha)/gamma(6) .* tt.^5 ./ T.^(5+alpha);
q = cap + ye.^2;
if nargin<3
    if isscalar(y)
        out = -y.^2 + q(1);
    else
        yy=y;
        if isrow(yy) && numel(yy)==numel(tt), yy=yy(:); end
        out = -yy.^2 + q;
    end
else
    out = -2*y;
end
end