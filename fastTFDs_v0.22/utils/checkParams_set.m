%-------------------------------------------------------------------------------
% Check and return arbitrary set from {k_1,k_2,...,k_J} samples
%
% could be for n (time) sample points also.
%
%-------------------------------------------------------------------------------
function k=checkParams_set(dec_param,N)
if(length(dec_param)>N)
  error('Frequency decimation parameters greater than N');
end
if(length(dec_param)==1)
  k=0:dec_param:N-1;
elseif(length(dec_param)>1)  
  k=dec_param;
else
  k=0:N-1;
end
if(max(k)>=N || min(k)<0) 
  error('Frequency decimation values need to within 0<=k<=N-1');
end


