function v10_benchmark_case
% Homogeneous MATLAB benchmark for official FLMM2/FHBVM/FHBVM2 and P2/P3/P4.
% All solvers are timed in the same MATLAB process and precision. The
% manufactured solution is smooth and the Pouzet starts use the computable
% truncation y=1; no exact startup values are injected.
more off; warning('off','all'); format long g;
try, maxNumCompThreads(1); catch, end %#ok<MAXNCT>
global V10_ALPHA V10_T
V10_ALPHA = str2double(getenv('ALPHA'));
V10_T = str2double(getenv('FINAL_TIME'));
if isnan(V10_ALPHA), V10_ALPHA=0.5; end
if isnan(V10_T), V10_T=1; end
reps = str2double(getenv('REPEATS')); if isnan(reps), reps=5; end
levels = [32 64 128 256];
outfile = sprintf('benchmark_a%.2f_T%g.csv',V10_ALPHA,V10_T);
fid=fopen(outfile,'w');
fprintf(fid,['alpha,T,method,level,parameter,nout,error,cpu_median,cpu_q1,' ...
             'cpu_q3,cpu_iqr,solver_reported_time,nfev,status,matlab_version\n']);
methods={'P2--PI2','P3--PI3','P4--PI4','FLMM2-BDF','FHBVM','FHBVM2'};
for im=1:numel(methods)
  method=methods{im};
  for ilev=1:numel(levels)
    level=levels(ilev);
    [t,y,nfev,status,times,reported,param]=run_repeated(method,level,ilev,reps);
    if isempty(t) || isempty(y) || any(~isfinite(y(:)))
      err=NaN; nout=0;
    else
      ye=v10_exact_smooth(t(:)); yy=y(:);
      if numel(yy)~=numel(ye)
        err=NaN; status=['fail:output-size-' num2str(numel(yy)) '-vs-' num2str(numel(ye))];
      else
        err=max(abs(yy-ye)./(1+abs(ye))); nout=numel(t);
      end
    end
    st=sort(times(isfinite(times)));
    if isempty(st)
      med=NaN; q1=NaN; q3=NaN; iq=NaN;
    else
      med=median(st); q1=local_quantile(st,.25); q3=local_quantile(st,.75); iq=q3-q1;
    end
    safe_status=strrep(strrep(status,',',';'),sprintf('\n'),' ');
    safe_version=strrep(version,',',';');
    fprintf(fid,['%.17g,%.17g,%s,%d,%.17g,%d,%.17g,%.17g,%.17g,%.17g,' ...
                 '%.17g,%.17g,%.17g,%s,%s\n'], ...
      V10_ALPHA,V10_T,method,level,param,nout,err,med,q1,q3,iq,reported,nfev,...
      safe_status,safe_version);
    fprintf('%s alpha=%.2f T=%g level=%d param=%g err=%.3e cpu=%.4g [%s]\n',...
      method,V10_ALPHA,V10_T,level,param,err,med,status);
  end
end
fclose(fid);
end

function q=local_quantile(x,p)
x=sort(x(:)); n=numel(x); pos=1+(n-1)*p; lo=floor(pos); hi=ceil(pos);
if lo==hi, q=x(lo); else, q=x(lo)+(pos-lo)*(x(hi)-x(lo)); end
end

function [t,y,nfev,status,times,reported,param]=run_repeated(method,level,ilev,reps)
try
  [t,y,nfev,reported,param]=run_method(method,level,ilev); status='ok';
catch err
  t=[]; y=[]; nfev=NaN; reported=NaN; param=NaN;
  status=['fail:' regexprep(err.message,',',';')];
end
times=NaN(1,reps);
if strcmp(status,'ok')
  % One untimed warm-up has already been executed above.
  repvals=NaN(1,reps);
  for k=1:reps
    try
      tic; [t,y,nfev,rt,param]=run_method(method,level,ilev); times(k)=toc; repvals(k)=rt;
    catch err
      status=['fail:' regexprep(err.message,',',';')]; break
    end
  end
  rr=repvals(isfinite(repvals)); if ~isempty(rr), reported=median(rr); end
end
end

function [t,y,nfev,reported,param]=run_method(method,level,ilev)
global V10_ALPHA V10_T
reported=NaN;
switch method
 case 'P2--PI2'
  param=level; [t,y,nfev]=pouzet_uniform(2,param);
 case 'P3--PI3'
  param=level; [t,y,nfev]=pouzet_uniform(3,param);
 case 'P4--PI4'
  param=level; [t,y,nfev]=pouzet_uniform(4,param);
 case 'FLMM2-BDF'
  param=level; h=V10_T/param;
  [t0,y0]=flmm2(V10_ALPHA,@(tt,yy)v10_problem_smooth(tt,yy),...
      @(tt,yy)v10_problem_smooth(tt,yy,1),0,V10_T,1,h,[],3,1e-12,100);
  t=t0(:); y=y0(:); nfev=NaN;
 case 'FHBVM'
  Mvals=[3 4 6 8]; param=Mvals(ilev);
  [t0,y0,stats]=fhbvm(@v10_problem_smooth,1,V10_T,param);
  t=t0(:); y=y0(:); nfev=NaN; reported=sum(stats(1:min(2,numel(stats))));
 case 'FHBVM2'
  param=level;
  [t0,y0,etim]=fhbvm2(@v10_problem_smooth,1,V10_T,param,1,1);
  t=t0(:); y=y0(:); nfev=NaN; reported=etim;
 otherwise
  error('unknown method')
end
end

function [t,y,nfev]=pouzet_uniform(stages,N)
global V10_ALPHA V10_T
alpha=V10_ALPHA; T=V10_T; q=stages; h=T/N;
[c,A,b]=tableau(alpha,stages);
t=linspace(0,T,N+1).'; y=zeros(N+1,1); rv=zeros(N+1,1);
% Computable startup: exact solution is 1+O(t^(5+alpha)); y=1 has an
% error O(h^(5+alpha)), above all target orders 1+s*alpha, s<=4.
y(1:min(q+1,N+1))=1;
for j=1:min(q+1,N+1), rv(j)=v10_problem_smooth(t(j),y(j)); end
nfev=min(q+1,N+1);
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
  off=idx-p; Vinv=local_vinv(off);
  coef=Vinv*vals(idx+1);
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

function Vinv=local_vinv(off)
persistent keys vals
if isempty(keys), keys={}; vals={}; end
key=sprintf('%g,',off);
for j=1:numel(keys)
  if strcmp(keys{j},key), Vinv=vals{j}; return; end
end
q=numel(off)-1; V=zeros(q+1,q+1);
for ii=1:q+1, V(ii,:)=off(ii).^(0:q); end
Vinv=V\eye(q+1); keys{end+1}=key; vals{end+1}=Vinv;
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