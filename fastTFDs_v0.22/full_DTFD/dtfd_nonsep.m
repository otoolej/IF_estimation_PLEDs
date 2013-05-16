%--------------------------------------------------------------------------------
% TFD (time-frequency distribution)
% Kernel Type:  nonseparable
%
%
% USE: [tf,R,g]=dtfd_nonsep(z,tfd_type,tfd_params)
%
% INPUT: 
%      z = input signal (either real-valued signal of length-N or
%          complex-valued analytic signal of length-2N)
%
%      tfd_type = { 'wvd' | 'swvd' | 'pwvd' | 'sep' | 'cw' | 'mb' }
%            wvd  - Wigner-Ville distribution
%            swvd - Smoothed Wigner-Ville distribution
%                   (lag-independent kernel)
%            pwvd - Pseudo Wigner-Ville distribution
%                   (Doppler-independent kernel)
%            sep  - Separable-kernel distribution 
%                   (combintation of SWVD and PWVD)
%            cw   - Choi-Williams distribution
%            mb   - Modified-B distribution
% 
%      tfd_params = cell of kernel parameters:
%            wvd  - {}
%            swvd - {win_length,win_type,[win_param]}
%                   e.g. {11,'hamm'}
%            pwvd - {win_length,win_type,[win_param]}
%                   e.g. {200,'cosh',0.1
%            sep  - { {win1_length,win1_type,[win1_param]}, 
%                    {win2_length,win2_type,[win2_param]}
%                   where win1 is the doppler window and win2 is the 
%                   lag window, e.g. { {11,'hamm'}, {200,'cosh',0.1} }
%            cw   - {sigma_parameter}
%            mb   - {beta_parameter} in the range 1<beta<0
%
% OUTPUT: 
%      tf = 2N x N TFD
%      R  = time-lag signal function
%      g  = kernel
%
% EXAMPLE:
%      N=128; 
%      x=gen_LFM(N,0.1,0.3)+gen_LFM(N,0.4,0.1);
%
%      c=dtfd_nonsep(x,'cw',{100}); 
%      vtfd(c,x);

%
%  Copyright (c) 2010, John M. O' Toole, The University of Queensland
%  All rights reserved.
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
%
%--------------------------------------------------------------------------------
function [tf,R,g]=dtfd_nonsep(z,tfd_type,tfd_params)
[z,N2,N,Nh]=check_analytic_signal(z);
tf=[]; R=[]; g=[];

global_vars;
set_DBvars('all',0); set_DBvars('DBtfd',0);




%-------------------------------------------------------------------------
% 0. Get Doppler-lag (smoothing) kernel
%-------------------------------------------------------------------------
g=gen_kernel(tfd_type,tfd_params,N);
g=g(:,1:(N+1));


%-------------------------------------------------------------------------
% 1. Get TL function K[n,m]
%    (K is the shitfed down matrix)
%-------------------------------------------------------------------------
K=zeros(N,N+1);
m=0:Nh-1; 
for n=0:N-1
  i1=mod(n+m,N2); i2=mod(n-m,N2); i3=mod(n+m+1,N2); 
  
  K(n+1,2.*m+1) = z(i1+1).*conj(z(i2+1));  % for m even
  K(n+1,2.*m+2) = z(i3+1).*conj(z(i2+1));  % for m odd.
end


%-------------------------------------------------------------------------
% 2. DFT to Doppler-lag domain
%-------------------------------------------------------------------------
af=fft_complex(K);


%-------------------------------------------------------------------------
% 3. Multiple kernel and signal function in Doppler-lag domain
%-------------------------------------------------------------------------
gaf=af.*g;


%-------------------------------------------------------------------------
% 4. DFT back to time-lag domain.
%-------------------------------------------------------------------------
R=ifft_complex(gaf);


%-------------------------------------------------------------------------
% 5.  Expand R for positive and negative lag values
%-------------------------------------------------------------------------
R=get_conjR(R,N2,N,Nh);


%-------------------------------------------------------------------------
% 6.  Go to the TF domain..
%-------------------------------------------------------------------------
tf=dft_Rdown(R,N2,N);
tf=tf./N2;  


tf=istfd_real(tf,z,mfilename);





%--------------------------------------------------------------------------------
% Expand time--lag function R[n,m]
%
%          from: 0<tau<(NT+1)
%            to: |tau|<NT
%
% R matrix has (T,T) discrete grid.
%--------------------------------------------------------------------------------
function Rnew=get_conjR(R,N2,N,Nh);
Nl=size(R,2); Nlh=ceil( Nl/2 );

Rnew=zeros(N,N2);
Rnew(:,1:Nl)=R(:,:);


% for even m ...
m=1:Nh-1;
Rnew(:,N2-2.*m+1)=conj( R(:,2.*m+1) );
% for odd m...
m=0:Nh-1;
Rnew(:,N2-2.*m)=conj( R(:,2.*m+2) );






