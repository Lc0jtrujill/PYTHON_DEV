function v10_benchmark_case
% Homogeneous GNU Octave benchmark for official FLMM2/FHBVM2 and P2/P3/P4.
more off; warning('off','all'); format long g;
global V10_ALPHA V10_T V10_LAMBDA V10_MU
V10_ALPHA = str2double(getenv('ALPHA'));
V10_T = str2double(getenv('FINAL_TIME'));
if isnan(V10_ALPHA), V10_ALPHA=0.5; end
if isnan(V10_T), V10_T=1; end
V10_LAMBDA=-0.2; V10_MU=0.1;
reps = str2double(getenv('REPEATS')); if isnan(reps), reps=5; end
Ns = [32 64 128 256];
outfile = sprintf('benchmark_a%.2f_T%g.csv',V10_ALPHA,V10_T);
fid=fopen(outfile,'w');
fprintf(fid,'alpha,T,method,N,nout,error,cpu_median,cpu_iqr,nfev,status,octave_version\n');
methods={'P2--PI2','P3--PI3','P4--PI4','FLMM2-BDF','FHBVM2'};
for im=1:numel(methods)
  method=methods{im};
  for N=Ns
    [t,y,nfev,status,times]=run_repeated(method,N,reps);
    if isempty(t) || any(~isfinite(y(:)))
      err=NaN; nout=0;
    else
      ye=v10_exact_smooth(t(:)); yy=y(:);
      err=max(abs(yy-ye)./(1+abs(ye))); nout=numel(t);
    end
    st=sort(times(isfinite(times)));
    if isempty(st), med=NaN; iq=NaN;
    else
      med=median(st); iq=st(max(1,ceil(.75*numel(st))))-st(max(1,ceil(.25*numel(st))));
    end
    fprintf(fid,'%.17g,%.17g,%s,%d,%d,%.17g,%.17g,%.17g,%.17g,%s,%s\n',...
      V10_ALPHA,V10_T,method,N,nout,err,med,iq,nfev,status,version);
    fprintf('%s alpha=%.2f T=%g N=%d err=%.3e cpu=%.4g [%s]\n',method,V10_ALPHA,V10_T,N,err,med,status);
  end
end
fclose(fid);
end

function [t,y,nfev,status,times]=run_repeated(method,N,reps)
try
  [t,y,nfev]=run_method(method,N); status='ok';
catch err
  t=[];y=[];nfev=NaN;status=['fail:' regexprep(err.message,',',';')];
end
times=NaN(1,reps);
if strcmp(status,'ok')
  for k=1:reps
    try
      tic; [t,y,nfev]=run_method(method,N); times(k)=toc;
    catch err
      status=['fail:' regexprep(err.message,',',';')]; break
    end
  end
end
end

function [t,y,nfev]=run_method(method,N)
global V10_ALPHA V10_T
switch method
 case 'P2--PI2'
  [t,y,nfev]=pouzet_uniform(2,N);
 case 'P3--PI3'
  [t,y,nfev]=pouzet_uniform(3,N);
 case 'P4--PI4'
  [t,y,nfev]=pouzet_uniform(4,N);
 case 'FLMM2-BDF'
  h=V10_T/N;
  [t0,y0]=flmm2(V10_ALPHA,@(tt,yy)v10_problem_smooth(tt,yy),...
      @(tt,yy)v10_problem_smooth(tt,yy,1),0,V10_T,1,h,[],3,1e-12,100);
  t=t0(:); y=y0(:); nfev=NaN;
 case 'FHBVM2'
  [t0,y0]=fhbvm2(@v10_problem_smooth,1,V10_T,N,1,1);
  t=t0(:); y=y0(:); nfev=NaN;
 otherwise
  error('unknown method')
end
end

function [t,y,nfev]=pouzet_uniform(stages,N)
global V10_ALPHA V10_T
alpha=V10_ALPHA; T=V10_T; q=stages; h=T/N;
[c,A,b]=tableau(alpha,stages);
t=linspace(0,T,N+1).'; y=zeros(N+1,1); rv=zeros(N+1,1);
y(1:q+1)=v10_exact_smooth(t(1:q+1));
for j=1:q+1, rv(j)=v10_problem_smooth(t(j),y(j)); end
nfev=q+1;
for n=q:N-1
  K=zeros(stages,1); K(1)=rv(n+1);
  for i=2:stages
    ti=t(n+1)+c(i)*h;
    H=1+pi_history(alpha,q,n,c(i),h,rv(1:n+1));
    Yi=H+h^alpha*(A(i,1:i-1)*K(1:i-1));
    K(i)=v10_problem_smooth(ti,Yi); nfev=nfev+1;
  end
  H=1+pi_history(alpha,q,n,1,h,rv(1:n+1));
  y(n+2)=H+h^alpha*(b*K);
  rv(n+2)=v10_problem_smooth(t(n+2),y(n+2)); nfev=nfev+1;
end
end

function H=pi_history(alpha,q,n,c,h,vals)
if n<=0, H=0; return; end
H=0;
for p=0:n-1
  if p<q, idx=0:q; else, idx=(p-q+1):(p+1); end
  off=idx-p; V=zeros(q+1,q+1);
  for ii=1:q+1, V(ii,:)=off(ii).^(0:q); end
  coef=V\vals(idx+1);
  AA=n+c-p; BB=AA-1; mom=zeros(q+1,1);
  for k=0:q
    mk=0;
    for ell=0:k
      mk=mk+nchoosek(k,ell)*AA^(k-ell)*(-1)^ell*...
        (AA^(alpha+ell)-BB^(alpha+ell))/(alpha+ell);
    end
    mom(k+1)=mk;
  end
  H=H+coef.'*mom;
end
H=h^alpha/gamma(alpha)*H;
end

function [c,A,b]=tableau(alpha,s)
if s==2
  c=[0;.5]; A=zeros(2); A(2,1)=c(2)^alpha/gamma(1+alpha);
  b2=1/(c(2)*gamma(2+alpha)); b=[1/gamma(1+alpha)-b2,b2];
elseif s==3
  c2=2/((1+alpha)*(2+alpha)); c3=2/(2+alpha); c=[0;c2;c3];
  A=zeros(3); A(2,1)=c2^alpha/gamma(1+alpha); A(3,2)=c3^alpha/gamma(1+alpha);
  b=[alpha/(2*(1+alpha)*gamma(1+alpha)),0,(2+alpha)/(2*(1+alpha)*gamma(1+alpha))];
else
  [c,A,b]=p4_tableau(alpha);
end
end

function [cv,A,b]=p4_tableau(alpha)
lo=2/(3+alpha); hi=2/(2+alpha); xx=linspace(lo+1e-12,hi-1e-12,4097);
ff=arrayfun(@(x)p4eq(alpha,x),xx); root=NaN;
for j=1:numel(xx)-1
  if ff(j)==0 || ff(j)*ff(j+1)<0
    root=fzero(@(x)p4eq(alpha,x),[xx(j),xx(j+1)]); break
  end
end
if isnan(root), error('P4 root not found'); end
c=root; d=p4d(alpha,c); mu1=1/gamma(2+alpha); mu2=2/gamma(3+alpha);
B=(mu1*d-mu2)/(c*(d-c)); b4=(mu2-mu1*c)/(d*(d-c)); b3=B/(1+alpha); b2=alpha*b3;
b1=1/gamma(1+alpha)-B-b4; a21=c^alpha/gamma(1+alpha); a32=a21;
S=d^(1+alpha)/(c*gamma(2+alpha)); a43=(S+(b3/b4)*a32)/(1+alpha); a42=S-a43; a41=d^alpha/gamma(1+alpha)-S;
cv=[0;c;c;d]; A=zeros(4); A(2,1)=a21; A(3,2)=a32; A(4,1:3)=[a41,a42,a43]; b=[b1,b2,b3,b4];
end
function d=p4d(a,c), d=(6/(3+a)-2*c)/(2-(2+a)*c); end
function z=p4eq(a,c)
d=p4d(a,c); z=a*c^(1+a)*((2+a)*d-2)+d^a*(2-(2+a)*c)*((2+a)*c-2*d);
end
