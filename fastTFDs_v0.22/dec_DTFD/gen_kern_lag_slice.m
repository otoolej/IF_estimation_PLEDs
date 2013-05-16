%-------------------------------------------------------------------------------
% Doppler-lag kernels for time--frequency distributions
% Return only a lag-slice of the kernel (see gen_kernel.m for full kernel)
%
%
%
% USE: g=gen_kern_lag_slice( lag_slice, kernel_type, kernel_params, N )
%
% INPUT: 
%   lag_slice     = lag value
%   kernel_type   = different kernels:
%                    'cw' - kernel for Choi-Williams distribution
%   kernel_params = cell of kernel parameters:
%                     cw - {sigma_parameter}
%   N             = signal length.
%
% OUTPUT:
%   g = Doppler--lag smoothing kernel, at one lag-value only
%
% EXAMPLE:
%      N=128;  lag_slice=3;
%      g=gen_kern_lag_slice( lag_slice, 'cw', {111}, 128 ); 
%
%      dopps=(-ceil(N/2)+1):ceil(N/2); 
%      clf; plot(dopps,shiftWin(g));
%      xlabel('Doppler');

% Copyright (C) 2010 John M. O' Toole, The University of Queensland
%
%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following
%  conditions are met:
%      * Redistributions of source code must retain the above
%        copyright notice, this list of conditions and the following
%        disclaimer.
%      * Redistributions in binary form must reproduce the above
%        copyright notice, this list of conditions and the following
%        disclaimer in the documentation and/or other materials
%        provided with the distribution.
%      * Neither the name of the The University of Queensland nor the 
%        names of its contributors may be used to endorse or promote 
%        products derived from this software without specific prior 
%        written permission.
%  
%  THIS SOFTWARE IS PROVIDED BY JOHN M. O' TOOLE ''AS IS'' AND ANY
%  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
%  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL JOHN M. O' TOOLE BE
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
%  OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
%  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
%  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
%  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
%  DAMAGE.
%-------------------------------------------------------------------------------
function g=gen_kern_lag_slice( lag_slice, kernel_type, kernel_params, N, G1 )
DB=0;
if(nargin<4) error('Need 3 input arguments.'); end
if(nargin<5 || isempty(G1)) G1=[]; end

N2=N*2;  

lag_sample_rate=1;
doppler_sample_rate=1/N;
g_lagslice=zeros(N,1);



g=get_kern(g_lagslice,lag_slice,kernel_type,kernel_params, ...
             doppler_sample_rate,lag_sample_rate,N2,N,G1);



% All kernels are real valued. 
g=real(g);



function g=get_kern(g,m,kernel_type,kernel_params,doppler_sample_rate, ...
                      lag_sample_rate,N2,N,G1)
Nd=size(g,1);
d_hlf=floor(Nd/2);
l=length(kernel_params);


switch kernel_type

  %---------------------------------------------------------------------
  % Wigner-Ville distribution.  For testing only
  %---------------------------------------------------------------------
 case {'wvd'}
  g=1;

  
  %---------------------------------------------------------------------
  % Choi-Williams
  %
  % g(l/NT,mT) = exp( -(2 pi m l/N)^2/ sigma ) 
  %
  %---------------------------------------------------------------------
 case {'cw', 'choi-williams', 'choi-will'}
  
  if(l>=1) 
    sigma=kernel_params{1};
  else
    sigma=11;
  end
  
  g(1) = 1;
  
  if(m==N)
      g=0;
  else
      if(m>N) m=2*N-m; end

      const = ((2*pi)^2)/sigma;
      u=1:d_hlf;
      u1=u.*doppler_sample_rate; m1=m*lag_sample_rate;
   
      g(u+1) = exp( -const.*(u1.*m1).^2 ).';
      g(Nd-u+1)=g(u+1); 
  end
  
  %---------------------------------------------------------------------
  % Product kernel (sometimes called a RID kernel)
  %
  % g(l/NT,mT)=H(lm/N)
  %---------------------------------------------------------------------
  case { 'prod', 'RID', 'product' }
    
    g(1) = 1;
  
    if(m==N)
        g=0;
    else
        if(m>N) m=2*N-m; end

        u=1:d_hlf;
        im=mod(u*m,length(G1));
        
        g(u+1)=G1(im+1);
        g(Nd-u+1)=g(u+1); 
    end
    

 otherwise
  error(['Unknown kernel type: ' kernel_type]);
  
end




