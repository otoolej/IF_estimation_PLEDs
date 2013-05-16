%-------------------------------------------------------------------------------
% Frequency decimation:
%
%    Y[k] = Y[k] + sum_{p=0}^{a-l} X[pJ+k] e^{j2c\pi(k+pJ)/N}
%
%    for k=0,1,...,J-1
%
% which returns IDFT{Y[k]} = y[an+c]
%
% where a and c are constants.
%
%-------------------------------------------------------------------------------
function y = fold_sig_full_odd( x, J, N, a, c )
% Assume that all parameters (J,a,c) are sane.

y=zeros(J,1);
x=x(:);


const=j*2*c*pi/N;

n=0:J-1;
y(n+1)=x(n+1).*exp(const.*n(:));   

for p=1:a-1
  ppi=p*J+n;
  y(n+1)=y(n+1)+(x(ppi+1).*exp(const.*ppi(:)));    
end

y=y./a;
