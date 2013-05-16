%-------------------------------------------------------------------------------
% Check and return L parameter;  a.n, for n=0,1,...,L-1
% where L=N/dec_param must be an integer
%
% could be for J (frequency) too.
%
%-------------------------------------------------------------------------------
function [L,dec_param]=checkParams_seq(dec_param,N)
if(length(dec_param)>1 || round(dec_param)~=dec_param)
  error(['dec_param parameter is integer scaler time decimation factor.']);
end
if(rem(N,dec_param)) 
  warning('2N (or N) not divisable by decimation factor. Increasing....');
  while( rem(N,dec_param) )  
    dec_param=dec_param+1;   
  end
  disp( ['... decimation factor now=' num2str(dec_param)]); 
end
L=N/dec_param;
if(L<=1) error('L (or J) is less than input signal length'); end
