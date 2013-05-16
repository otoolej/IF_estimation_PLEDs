%-------------------------------------------------------------------------------
% Frequency decimation:
%
%    y[n] = y[n] + sum_{p=0}^{a-l} x[pJ+n]
%
%    for n=0,1,...,Jh
%
%-------------------------------------------------------------------------------
function y = fold_sig_half( x, J, Jh, a )
% Assume that all parameters (J,Jh,a) are sane.

y=zeros(Jh+1,1);
x=x(:);

n=0:Jh;
for p=0:a-1
  y(n+1)=y(n+1)+x(p*J+n+1);    
end


