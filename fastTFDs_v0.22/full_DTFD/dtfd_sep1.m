%--------------------------------------------------------------------------------
% TFD (time-frequency distribution)
%
% Kernel Type: separable-kernel.  
% Using algorithm which starts in the time-lag domain, as follows:
%       time--lag -> Doppler--lag -> time--lag --> time--frequency
% See [1] for details.
%
% 
% USE: [tf,R,G1,g2]=dtfd_sep1(z,doppler_win_params,lag_win_params, ...
%                                 Ntime,Nfreq)
%
% INPUT: 
%      z = input signal (either real-valued signal of length-N or
%          complex-valued analytic signal of length-2N)
%
%      doppler_win_params = {win_length,win_type,[win_param],dopp_or_time} 
%                   with:
%                     win_length  - length (odd value) of window,
%                     win_type    - name of window, e.g. 'hamm' for
%                                   Hamming window
%                     win_param   - if window has a parameter value,
%                                   e.g. 0.2 for Tukey window
%                     dopp_or_time - either:
%                       0 = to define window in the Doppler-domain
%                       1 = to define window in the time-domain
%                   e.g. {11,'hamm',0,1} or {200,'cosh',0.1} 
%
%      lag_win_params = {win_length,win_type,[win_param],lag_of_freq} 
%                        e.g. {111,'hamm'} or {127,'tukey',0.2} 
%
%      Ntime = number of sample points in the time direction for 
%              the TFD, Ntime <= 2N
%      Nfreq = number of sample points in the frequency direction for
%              the TFD, Nfreq <= N
%      
%
% OUTPUT: 
%      tf = Ntime x N TFD
%      R  = time-lag signal function
%      G1 = lag-independent kernel 
%      g2 = Doppler-independent kernel 
%
% EXAMPLE:
%      N=10000; Ntime=256; Nfreq=256;
%      x=gen_LFM(N,0.1,0.3)+gen_LFM(N,0.4,0.1);
%
%      c=dtfd_sep1(x,{51,'hamm',0,1},{271,'hann'},Ntime,Nfreq); 
%      vtfd(c,x);

%   Copyright (c) 2010, John M. O' Toole, The University of Queensland
%   All rights reserved.
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
%--------------------------------------------------------------------------------
function [tf,Rdown,G1,g2]=dtfd_sep1(z,doppler_win_params,lag_win_params, ...
                             Ntime,Nfreq)
if(nargin<4) Ntime=[]; end
if(nargin<5) Nfreq=[]; end
[z,N2,N,Nh]=check_analytic_signal(z);
tf=[]; R=[]; g2=[]; G1=[];


global_vars; 
set_DBvars('all',0); set_DBvars('DBplot',0); set_DBvars('DBtfd',0);
DBprofile=0;


%-------------------------------------------------------------------------
% 1. Get kernel g(nu,tau) = G_1(nu)  (defined in the EITHER time domain 
%     or doppler domain)
%
%  extent: 0<nu<1/T;  G_1(nu) is length Q.
%-------------------------------------------------------------------------
[G1,Q,time_doppler,Ntime]=get_Doppler_kernel(doppler_win_params,N, ...
                                                          Ntime);
% If Ntime is odd then adjust and check out G1 again: (should
% really do this check in 'get_Doppler_kernel.m' function)
while(rem(Ntime,2)==1)
  Ntime=Ntime+1;
  [G1,Q,time_doppler,Ntime]=get_Doppler_kernel(doppler_win_params,N, ...
                                                          Ntime);
end
Qh=floor(Q/2);
Nh_time=floor(Ntime/2);
Q_even_factor=rem(Q,2)-1;
if(FAST_DTFD_DBvars) dispVars(Ntime,Nh_time,Q,Qh,Q_even_factor); end


% Modulate Doppler-window
if(Nh_time~=N)
  G1mod=zeros(Q,1);
  lp=0:Qh; ln=1:(Qh+Q_even_factor);
  G1mod(lp+1)=G1(lp+1).*exp(-j*pi/N.*lp(:)).*exp(j*pi/Nh_time.*lp(:));
  G1mod(Q-ln+1)=G1(Q-ln+1).*exp(-j*pi/N.*(N-ln(:))).* ...
      exp(j*pi/Nh_time.*(Nh_time-ln(:)));
else
  G1mod=G1;
end



%-------------------------------------------------------------------------
% 2. Get lag window g_2[m]
%-------------------------------------------------------------------------
[g2,P,Ph,Nfreq]=get_lag_kernel(lag_win_params,N2,N,Nfreq);
g2=g2.';
N2_freq=2*Nfreq;
Nh_freq=ceil(Nfreq/2);
Pq=ceil(Ph/2);
Ph_even_factor=rem(Ph,2)-1;
if(FAST_DTFD_DBvars) dispVars(Nfreq,Nh_freq,P,Ph,Pq,Ph_even_factor); end




R=zeros(Nh_time,Ph+1);

% Set some indicies:
n=0:N-1;
l=0:Nh_time-1;
lp=0:Qh; ln=1:(Qh+Q_even_factor);
lz=(Qh+1):(Nh_time-Qh-Q_even_factor-1);


%-------------------------------------------------------------------------
% 3. Form time-lag function K[n,m] multiplied by g2[m]
%    (matrix Rdown[n,m])
%-------------------------------------------------------------------------
% 3A. For m even:
K_meven_slice=zeros(N,1); 
saf_meven_slice=zeros(Nh_time,1);    
for m=0:Pq-1-Ph_even_factor
  i1=mod(n+m,N2); i2=mod(n-m,N2);

  
  %-------------------------------------------------------------------------
  % a). Form time-lag function and window in the lag direction 
  %-------------------------------------------------------------------------
  K_meven_slice(n+1) = z(i1+1).*conj(z(i2+1)).*g2(2*m+1);   

  
  %-------------------------------------------------------------------------
  % b) Transform to Doppler-lag [l,m] domain
  %-------------------------------------------------------------------------
  af_meven_slice=fft_complex(K_meven_slice);


  %-------------------------------------------------------------------------
  % c). Multiply by G1[l]
  %-------------------------------------------------------------------------
  saf_meven_slice(lp+1)=af_meven_slice(lp+1).*G1(lp+1);
  saf_meven_slice(Nh_time-ln+1)=af_meven_slice(N-ln+1).*G1(Q-ln+1);
  saf_meven_slice(lz+1)=0;
  
  
  %-------------------------------------------------------------------------
  % d). Transform back to time-lag domain [n,m]
  %-------------------------------------------------------------------------
  R(:,2*m+1)=ifft_complex(saf_meven_slice);

  if(DBprofile); dispVars(m); end
end
clear K_meven_slice, saf_meven_slice;



% 3B. For m odd:
K_modd_slice=zeros(N,1); 
saf_modd_slice=zeros(Nh_time,1);    
for m=0:Pq-1
  i2=mod(n-m,N2); i3=mod(n+m+1,N2); 
  
  %-------------------------------------------------------------------------
  % a). Form time-lag function and window in the lag direction 
  %-------------------------------------------------------------------------
  K_modd_slice(n+1) = z(i3+1).*conj(z(i2+1)).*g2(2.*m+2);  
  

  %-------------------------------------------------------------------------
  % b) Transform to Doppler-lag (nu,tau) domain
  %
  %  extent: |nu|<1/2T and 0<tau<(NT+1)
  %-------------------------------------------------------------------------
  af_modd_slice=fft_complex(K_modd_slice);

  
  %-------------------------------------------------------------------------
  % c). Multiply by G1(l)
  %-------------------------------------------------------------------------
  saf_modd_slice(lp+1)=af_modd_slice(lp+1).*G1mod(lp+1);
  saf_modd_slice(Nh_time-ln+1)=af_modd_slice(N-ln+1).*G1mod(Q-ln+1);
  saf_modd_slice(lz+1)=0;
  

  %-------------------------------------------------------------------------
  % d). Transform back to time-lag domain (t,tau)
  %-------------------------------------------------------------------------
  R(:,2*m+2)=ifft_complex(saf_modd_slice);
  
  if(DBprofile); dispVars(m); end
end
clear K_modd_slice saf_modd_slice;



%-------------------------------------------------------------------------
% 5. Expand in the lag direction
%
%  extent: |t|<NT/2 and |tau|<NT/2
%-------------------------------------------------------------------------
if(DBprofile) dispVars('rearrange R matrix.'); end
m1=1:Nh_freq-1; m2=0:Nh_freq-1;
Rdown=zeros(Nh_time,N2_freq);
Rdown(:,1:(Ph+1))=R(:,1:(Ph+1));
Rdown(:,N2_freq-2.*m1+1)=conj(Rdown(:,2.*m1+1));
Rdown(:,N2_freq-2.*m2)=conj(Rdown(:,2.*m2+2));




%-------------------------------------------------------------------------
% 6. Finally, transform to the time-frequency domain
%-------------------------------------------------------------------------
if(DBprofile); dispVars('DFT to TF domain.'); end
tf=dft_Rdown(Rdown,Ntime,Nfreq);

%NOT sure on the scaling factor here??
% $$$ tf=scale_dtfd(tf,Nfreq);
tf=scale_dtfd(tf,Ntime);



tf=istfd_real(tf,z,mfilename);


