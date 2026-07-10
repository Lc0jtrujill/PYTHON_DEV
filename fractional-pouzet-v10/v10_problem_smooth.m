function out = v10_problem_smooth(t,y,flag)
% Chain-active manufactured scalar Caputo problem.
% Exact solution: y(t)=1+(t/T)^(1+alpha).
% The exact density is Gamma(2+alpha)*t/T^(1+alpha).
global V10_ALPHA V10_T V10_LAMBDA V10_MU
if nargin==0
    out = V10_ALPHA;
    return
end
alpha = V10_ALPHA; T = V10_T;
tt = t(:);
ye = 1 + (tt./T).^(1+alpha);
if nargin<3
    g = gamma(2+alpha).*tt./T.^(1+alpha);
    if isscalar(y)
        e = y-ye(1);
        out = g(1) + V10_LAMBDA.*e + V10_MU.*e.^2;
    else
        yy = y;
        if isrow(yy) && numel(yy)==numel(tt), yy = yy(:); end
        e = yy-ye;
        out = g + V10_LAMBDA.*e + V10_MU.*e.^2;
    end
else
    if isscalar(y)
        e = y-ye(1);
        out = V10_LAMBDA + 2*V10_MU.*e;
    else
        yy = y;
        if isrow(yy) && numel(yy)==numel(tt), yy = yy(:); end
        e = yy-ye;
        out = V10_LAMBDA + 2*V10_MU.*e;
    end
end
end
