%-------------------------------------------------------------------------------
% Doppler-lag kernels for time--frequency distributions
%
% USE: g=gen_kernel( kernel_type, kernel_params, N )
%
% INPUT: 
%
%   kernel_type = { 'wvd' | 'swvd' | 'pwvd' | 'sep' | 'cw' | 'mb'}
%            wvd  - kernel for Wigner-Ville distribution
%            swvd - kernel for smoothed Wigner-Ville distribution
%                   (lag-independent kernel)
%            pwvd - kernel for pseudo Wigner-Ville distribution
%                   (Doppler-independent kernel)
%            sep  - kernel for separable kernel (combintation of SWVD and PWVD)
%            cw   - kernel for Choi-Williams distribution
%            mb   - kernel for Modified-B distribution
% 
%   kernel_params = cell of kernel parameters:
%            wvd  - {}
%            swvd - {win_length,win_type,[win_param],[Doppler_or_time]}
%                   e.g. {11,'hamm',0,1}
%                   the parameter Doppler_or_time is either
%                       0 = to define window in the time-domain
%                       1 = to define window in the Doppler-domain
%            pwvd - {win_length,win_type,[win_param]}
%                   e.g. {200,'cosh',0.1}
%            sep  - { {win1_length,win1_type,[win1_param]}, 
%                    {win2_length,win2_type,[win2_param]}
%                   where win1 is the Doppler window and win2 is the lag window.
%                   e.g. { {11,'hamm'}, {200,'cosh',0.1} }
%            cw   - {sigma_parameter}
%            mb   - {beta_parameter} in the range 1<beta<0
%
%   N = signal length.
%
% OUTPUT:
%   g = Doppler--lag smoothing kernel
%
%
% EXAMPLE
%      N=128; 
%      g=gen_kernel( 'sep',{{N-1,'cosh',0.1,1},{51,'hamm'}},N); 
%      clf; mesh( fftshift(g) );
%      xlabel('lag'); ylabel('Doppler');
%      title('Separable kernel');
%


% Copyright (C) 2007,2008 John M. O' Toole, The University of Queensland
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
function g=gen_kernel( kernel_type, kernel_params, N)
if(nargin<3) error('Need 3 input arguments.'); end
N2=N*2;  

lag_sample=1;
doppler_sample=1/N;
g=zeros(N,N2);


g=get_kern(g,kernel_type,kernel_params,doppler_sample,lag_sample,N2,N);


% All kernels are real valued. 
g=real(g);





function g=get_kern(g,kernel_type,kernel_params,doppler_sample,lag_sample,N2,N)
%-------------------------------------------------------------------------------
% Generate Doppler--lag kernel function
%-------------------------------------------------------------------------------
[Nd,Nl]=size(g);
d_hlf=floor(Nd/2);
l_hlf=floor((Nl-1)/2);
l=length(kernel_params);


switch kernel_type
  
 case {'wvd'}
  g(:,:)=1;

  
  %---------------------------------------------------------------------
  % Smoothed Wigner-Ville (Lag Independent (LI) kernel)
  % g(l/NT,mT) = W(l/NT)
  %---------------------------------------------------------------------
  case 'swvd'
    G1=get_Doppler_kernel(kernel_params,N,2*N);
    if(length(kernel_params)>3 & kernel_params{4}==1)
        G1=padWin(G1,N);
    end
  
    for m=1:Nl
        g(:,m)=G1;
    end
  
  
  %---------------------------------------------------------------------
  % Pseudo-Wigner-Ville (Doppler Independent (DI) kernel)
  % g(l/NT,mT) = g2(mT)
  %---------------------------------------------------------------------
 case 'pwvd'
  g2=get_lag_kernel(kernel_params,2*N,N,Nl);
  g2=padWin(g2,Nl);
   
  for l=0:Nd-1
      g(l+1,:)=g2;
  end

  %---------------------------------------------------------------------
  % Seperable Kernel
  %
  % g(l/NT,mT) = G1(l/NT)g2(mT) 
  %  
  %---------------------------------------------------------------------
 case { 'sep', 'sep-full' }

  if(l<2)
    error('Need at least two windows parameters.'); 
  end
    
  g1=get_kern(g,'swvd',kernel_params{1},doppler_sample,lag_sample,N2,N);
  g2=get_kern(g,'pwvd',kernel_params{2},doppler_sample,lag_sample,N2,N);  

  g=g1.*g2;

  
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
  
  
  g(1,1:Nl) = 1;
  g(1:Nd,1) = 1;
  g(1,l_hlf+2) = 0;
  
   const = ((2*pi)^2)/sigma;
   u=1:d_hlf;
   for m=1:l_hlf
     u1=u.*doppler_sample; m1=m*lag_sample;
     g(u+1,m+1) = exp( -const.*(u1.*m1).^2 ).';
     
     g(Nd - u + 1, m + 1) = g(u + 1, m + 1); 
     g(u + 1, Nl - m + 1) = g(u + 1, m + 1); 
     g(Nd - u + 1, Nl - m + 1) = g(u + 1, m + 1); 
   end
   
  %---------------------------------------------------------------------
  % Modified-B (a specific SWVD or LAG INDEPENDENT kernel)
  %
  % G(nT,mT) = cosh^{-2*beta}(n)   (defined in time--lag domain)
  %---------------------------------------------------------------------
 case { 'mb', 'modified-b' }
  
  g=get_kern(g,'swvd',{N,'cosh',kernel_params{1}},doppler_sample,lag_sample,N2,N);


  %---------------------------------------------------------------------
  % Product kernel (sometimes called a RID kernel)
  %
  % G(nT,mT) = 1/|mT| h(n/|m|) or g(l/NT,mT)=H(lm/N)
  %---------------------------------------------------------------------
  case { 'prod', 'RID', 'product' }
    
    % oversample the window:
    L_dopp=Nd*Nl;
    
    G1=get_Doppler_kernel(kernel_params,L_dopp,2*N);
    G1=real(G1);
    
    g(1,1:Nl)=1;
    g(1:Nd,1)=1;
    g(1,l_hlf+2)=0;

    u=1:d_hlf;
    for m=1:l_hlf
        im=mod(u*m,L_dopp);
        
        g(u+1,m+1) = G1(im+1);
     
        g(Nd - u + 1, m + 1) = g(u + 1, m + 1); 
        g(u + 1, Nl - m + 1) = g(u + 1, m + 1); 
        g(Nd - u + 1, Nl - m + 1) = g(u + 1, m + 1); 
    end
  
 otherwise
  error(['Unknown kernel type: ' kernel_type]);
  
end


