%-------------------------------------------------------------------------------
% Decimated WVD (Wigner-Ville distribution).  See [1] for details.
%
%
% USE: wvd=dwvd_grid1(z,time_dec,freq_dec)
%
% INPUT:
%      z = N-point real-valued time-domain signal
%      time_dec  = constant where L=(2N/time_dec) is an integer.
%      freq_dec  = either a constant, where V=(N/freq_dec) is an integer
%                  OR
%                  set freq_dec=[k1,k2,...,kV], where v<N and 0<=ki<=N  
%
% OUTPUT:
%     wvd = L x V WVD, where
%                  - L=(2N/time_dec) 
%                  - V is the length of the set freq_dec (or if freq_dec is a scalar
%                    then V=N/freq_dec)
%     K   = time-lag signal function
%
% EXAMPLE:
%      N=1024; time_dec=8; 
%      freq_dec=[50:10:300,301:2:700,701:10:1000];
%      x=gen_LFM(N,0.1,0.3)+gen_LFM(N,0.4,0.1);
%
%      wc=dwvd_grid1(x,time_dec,freq_dec); 
%      clf; vtfd(wc);
% 


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
%-------------------------------------------------------------------------------
function [w,K] = dwvd_grid1(z,time_dec,freq_dec)
if(nargin<2) time_dec=1; end
if(nargin<3) freq_dec=1; end

[z,N2,N,Nh]=check_analytic_signal(z);
global_vars;
set_DBvars('all',0);  set_DBvars('dbtfd',0);
w=[]; K=[];



%-------------------------------------------------------------------------------
% 0. Check decimation parameters
%-------------------------------------------------------------------------------
k=checkParams_set(freq_dec,N);
[L,time_dec]=checkParams_seq(time_dec,N2);
if(FAST_DTFD_DBvars) dispVars(time_dec,freq_dec,L); end


%-------------------------------------------------------------------------------
% 1. Get Doppler-frequency function for K(:,kfreq) folded in the l direction for
% time-decimation.
%-------------------------------------------------------------------------------
K=get_Ksiaf(z,N2,N,Nh,k,L,time_dec); 


%-------------------------------------------------------------------------------
% 3. FFT to the TF domain.
%-------------------------------------------------------------------------------
w=ifft_conj_symm(K);
w=w./N2;

w=istfd_real(w,z,'decFreq_DWVD_C');



 



function K=get_Ksiaf(z,N2,N,Nh,k,L,time_dec)
%-------------------------------------------------------------------------
% 3. Multiple DF function K(l/2NT,k/2NT) by kernel G_1(l/2NT) 
%
%  extent: 0<nu<1/T;  g_1(t) is length Q.
%-------------------------------------------------------------------------
Z=fft(z); 
J=length(k);
Lh=ceil(L/2);
L_odd_factor=rem(L,2);

K=zeros(L,J);
Kslice=zeros(N2,1);

l=0:N;  li=0:Lh; lb=1:N-1; lbL=1:(Lh-1);
for ki=1:length(k)
  i1=mod(k(ki)+l,N2); i2=mod(k(ki)-l,N2); 
  
  Kslice(l+1) = Z(i1+1).*conj(Z(i2+1));  
  Kslice(N2-lb+1)=conj( Kslice(lb+1) ); 

  % If need to fold for time decimation then do so here:
  if(L~=N2)
    K(li+1,ki)=fold_sig_half(Kslice,L,Lh,time_dec);
    K(li+1,ki)=K(li+1,ki)./time_dec;
    K(L-lbL+1,ki)=conj( K(lbL+1,ki) );          
  else
    K(:,ki)=Kslice;
  end
end

