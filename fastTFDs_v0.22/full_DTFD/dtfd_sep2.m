%-----------------------------------------------------------------------------------------
% TFD (time-frequency distribution)
%
% Kernel Type: separable-kernel.  
% Using algorithm which starts in the time-lag domain, as follows:
%       Doppler--frequency -> Doppler--lag -> Doppler--frequency --> time--frequency
% See [1] for details.
%
% 
% USE: [tf,R,G1,g2]=dtfd_sep2(z,doppler_win_params,lag_win_params, ...
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
function [tf,R,G1,g2] = dtfd_sep2(z,doppler_win_params, ...
                                  lag_win_params,Ntime,Nfreq)
if(nargin<4) Ntime=[]; end
if(nargin<5) Nfreq=[]; end
[z,N2,N,Nh]=check_analytic_signal(z);
tf=[]; R=[]; g2=[]; G1=[];


global_vars; 
set_DBvars('all',0); set_DBvars('DBplot',0); set_DBvars('DBtfd',0);



%-------------------------------------------------------------------------
% 1. Get Doppler window G_1[l]  
%    (defined in the EITHER time domain or Doppler domain)
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
Q_even=rem(Q,2)-1;
if(FAST_DTFD_DBvars) dispVars(Ntime,Nh_time,Q,Qh,Q_even); end
if(FAST_DTFD_DBplot) plotc(G1); end
g=G1;

%-------------------------------------------------------------------------
% 2. Get lag window g_2[m]
%-------------------------------------------------------------------------
[g2,P,Ph,Nfreq]=get_lag_kernel(lag_win_params,N2,N,Nfreq);
g2=g2.';
N2_freq=2*Nfreq;
P_even=rem(P,2)-1;
if(FAST_DTFD_DBvars) dispVars(P,Ph,Nfreq); end
if(FAST_DTFD_DBplot) plot(g2); end

%-------------------------------------------------------------------------
% 3. Doppler-slice at a time
%-------------------------------------------------------------------------
R=zeros(Qh+1,N2_freq);
Rslice=zeros(N2,1);
saf_slice=zeros(N2_freq,1);
k=0:N2-1; 
m1=0:Ph; m2=1:Ph+P_even;

Z=fft(z);
for l=0:Qh
  %-------------------------------------------------------------------------
  % a) Doppler slice of Doppler-frequency function and window in
  %    frequency direction
  %-------------------------------------------------------------------------
  i1=mod(k+l,N2); i2=mod(k-l,N2); 
  Rslice(k+1) = Z(i1+1).*conj(Z(i2+1)).*G1(l+1);  
  
  %-------------------------------------------------------------------------
  % b) Transform to Doppler-lag domain
  %-------------------------------------------------------------------------
  af_slice=ifft_complex(Rslice);
  
  %-------------------------------------------------------------------------
  % c) Multiply by g2[m] 
  %-------------------------------------------------------------------------
  saf_slice(m1+1)=af_slice(m1+1).*g2(m1+1).';
  saf_slice(N2_freq-m2+1)=af_slice(N2-m2+1).*g2(P-m2+1).';  

  %-------------------------------------------------------------------------
  % d) Transform  back to Doppler-frequency domain
  %-------------------------------------------------------------------------
  R(l+1,:)=fft_complex(saf_slice);
end




%-------------------------------------------------------------------------
% 4a. Transform from K_across(l/2NT,k/2NT) 
%      from
%          0<nu<1/2T, |f|<1/2T
%      to 
%          0<nu<1/T, 0<f<1/2T
%-------------------------------------------------------------------------
Racross=K_move(R,Nh_time,Nfreq,Qh,Q_even);

%-------------------------------------------------------------------------
% 4b. Expand from in Doppler direction from
%          0<nu<1/T
%      to 
%          |nu|<1/T
%-------------------------------------------------------------------------
Racross=get_conj_dfR_long(Racross,Ntime,Nfreq,Nh_time);



%-------------------------------------------------------------------------
% 5. Go to time--frequency domain
%-------------------------------------------------------------------------
tf=ifft_conj_symm(Racross)./Ntime;


tf=istfd_real(tf,z,mfilename);







function R=get_RDLacross_wide(z,N2,G1,Qh)
%-------------------------------------------------------------------------
% 1. Get K_across(l/NT,k/2NT) , for l=0,1,..,Q/2-1; k=0,1,...,2N-1
%    
%    Extent: 0<nu<1/2T, |f|<1/2T
%-------------------------------------------------------------------------
Z=fft(z);
R=zeros(Qh+1,N2);
k=0:N2-1; 
for l=0:Qh
  i1=mod(k+l,N2); i2=mod(k-l,N2); 
  
  R(l+1,k+1) = Z(i1+1).*conj(Z(i2+1)).*G1(l+1);  
end


function R=get_RDLacross_wide_slice(Z,l,N2,G1,Qh)
%-------------------------------------------------------------------------
% 1. Get K_across(l/NT,k/2NT) , for k=0,1,...,2N-1
%-------------------------------------------------------------------------
Rslice=zeros(N2,1);
k=0:N2-1; 
for l=0:Qh
  i1=mod(k+l,N2); i2=mod(k-l,N2); 
  
  Rslice(l+1,k+1) = Z(i1+1).*conj(Z(i2+1)).*G1(l+1);  
end






function Kout=K_move(Kin,Nh_time,Nfreq,Qh,Q_even)
%--------------------------------------------------------------------------------
%
% Convert from (Ntime/4 x 2Nfreq) K_across(l/2NT,k/2NT) matrix 
% to the (Ntime/2 x Nfreq) K_across(l/2NT,k/2NT) one.
%
% That is, convert from 
%   0<nu<1/2T, |f|<1/2T 
% to 
%   0<nu<1/T,  0<f<1/2T 
%
%--------------------------------------------------------------------------------
Kout=zeros(Nh_time+1,Nfreq);

% Fill in the bits that share the same area
l=0:Qh;
for k=0:Nfreq-1
  Kout(l+1,k+1)=Kin(l+1,k+1);
end

% and then use the symmetry to patch in the other
% bits.
l=0:Qh+Q_even;
for k=0:Nfreq-1
  Kout(Nh_time-l+1,k+1)=conj( Kin(l+1,Nfreq+k+1) );
end





function Rnew=get_conj_dfR_long(R,Ntime,Nfreq,Nh_time);
%--------------------------------------------------------------------------------
% Expand Doppler--frequency function R[l,k]
%
%          from: 0<nu<1/T
%            to: |nu|<1/T
%
% R matrix has uniform (l/NT,k/2NT) discrete grid.
%--------------------------------------------------------------------------------
Rnew=zeros(Ntime,Nfreq);
Rnew(1:Nh_time+1,:)=R(:,:);
Ntime_even=rem(Ntime,2)-1;

global_vars;
if(FAST_DTFD_DBvars) dispVars(Ntime_even); end



l=1:Nh_time+Ntime_even;
Rnew(Ntime-l+1,:)=conj( R(l+1,:) );






