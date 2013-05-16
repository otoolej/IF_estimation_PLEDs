%-------------------------------------------------------------------------------
% Check input signal and return analytic associate if necessary
%
%
%-------------------------------------------------------------------------------
function [z,N2,N,Nh]=check_analytic_signal(z)
if(isreal(z)) 
  z=get_analytic(z); 
end
N2=length(z); 
N=N2/2; 
Nh=ceil(N/2);

if(N~=fix(N))
  error('Analytic signal must be length 2N.');
end
