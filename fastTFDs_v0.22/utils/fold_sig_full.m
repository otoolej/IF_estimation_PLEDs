%-------------------------------------------------------------------------------
% Frequency decimation:
%
%    y[n] = y[n] + sum_{p=0}^{a-l} x[pJ+n]
%
%    for n=0,1,...,J-1
%
% which returns DFT{y[n]} = Y[ak]
%
%-------------------------------------------------------------------------------
function y = fold_sig_full( x, J, a )
% Assume that all parameters (J,a) are sane.

y=zeros(J,1);
x=x(:);

n=0:J-1;
y(n+1)=x(n+1);    
for p=1:a-1
  y(n+1)=y(n+1)+x(p*J+n+1);    
end

% scale by decimation factor:
y=y./a;
