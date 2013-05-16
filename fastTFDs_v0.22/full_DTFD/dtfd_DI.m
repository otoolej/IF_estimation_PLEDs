%-------------------------------------------------------------------------------
% TFD (time-frequency distribution)
%
% Kernel Type: Doppler-independent kernel.  See [1] for details.
%
%
%  USE: [tf,R,g2] = dtfd_DI(z,tfd_params,Nfreq)
%
%  INPUT:
%      z          = input signal (either real-valued signal of length-N or
%                   complex-valued analytic signal of length-2N)
%      tfd_params = cell of {win_length,win_type,[win_param],[lag_or_freq]}
%                   with:
%                     win_length  - length (odd value) of window,
%                     win_type    - name of window, e.g. 'hamm' for
%                                   Hamming window
%                     win_param   - if window has a parameter value,
%                                   e.g. 0.2 for Tukey window
%                     lag_or_freq - either:
%                       0 = to define window in the lag-domain
%                       1 = to define window in the frequency-domain
%      Nfreq      = number of sample points in the frequency direction for
%                   the TFD, Nfreq <= N
%
%  OUTPUT:
%     tf = 2N x Nfreq TFD
%     R  = time-lag signal function
%     g2 = lag-independent kernel
%
%  EXAMPLE:
%     N=512; Nfreq=128;
%     x=gen_LFM(N,0.1,0.3)+gen_LFM(N,0.4,0.1);
%
%     tf=dtfd_DI(x,{151,'hamm'},Nfreq); 
%     clf; vtfd(tf,x);
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
function [tf,R,g2] = dtfd_DI(z,tfd_params,Nfreq)
if(nargin<3) Nfreq=[]; end
[z,N2,N,Nh]=check_analytic_signal(z);
tf=[]; R=[]; g2=[];


global_vars;
set_DBvars('all',0); set_DBvars('dbtfd',0);

%-------------------------------------------------------------------------
% 1. Get lag-window g_2[m] of length P
%-------------------------------------------------------------------------
set_DBvars('dbplot',0);
[g2,P,Ph,Nfreq]=get_lag_kernel(tfd_params,N2,N,Nfreq);
g2=g2(:);
Nh_freq=ceil(Nfreq/2);
Pq=ceil(Ph/2);

if(FAST_DTFD_DBvars) dispVars(P,Ph,Pq,Nfreq,Nh_freq); end
if(FAST_DTFD_DBplot) figure(334); plot(g2.'); end


%-------------------------------------------------------------------------------
% 2. Get smoothed TL function for K[n,m]g2[m], 
%    (K is the shitfed across matrix)
%-------------------------------------------------------------------------------
% First sort out even--odd time indices
R=zeros(N2,Nfreq);
Ph_even_factor=rem(Ph,2)-1;

%------------------------------------
% for n even values:
%------------------------------------
m=0:(Pq-1-Ph_even_factor); mh=1:Nh_freq-1;
for n=0:N-1
  i1=mod(n+m,N2); i2=mod(n-m,N2);
  
  R(2*n+1,m+1)=z(i1+1).*conj(z(i2+1)).*g2(2.*m+1); 
  R(2*n+1,Nfreq-mh+1)=conj( R(2*n+1,mh+1));  
end

%------------------------------------
% for n odd values:
%------------------------------------
m=0:(Pq-1); mh=0:Nh_freq-1;
for n=0:N-1  
  i3=mod(n+m+1,N2); i2=mod(n-m,N2);

  R(2*n+2,m+1)=z(i3+1).*conj(z(i2+1)).*g2(2.*m+2); 
  R(2*n+2,Nfreq-mh)=conj(R(2*n+2,mh+1));    
end



%-------------------------------------------------------------------------------
% 3. FFT to the TF domain.
%-------------------------------------------------------------------------------
tf=dft_Racross(R,N2,Nfreq,N);
tf=tf./N2;  

tf=istfd_real(tf,z,mfilename);



