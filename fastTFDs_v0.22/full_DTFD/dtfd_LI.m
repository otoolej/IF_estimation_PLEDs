%--------------------------------------------------------------------------------
% TFD (time-frequency distribution)
%
% Kernel Type: Doppler-independent kernel.  See [1] for details.
%
%
%  USE: [tf,R,G1] = dtfd_LI(z,tfd_params,Ntime)
%
%  INPUT:
%      z          = input signal (either real-valued signal of length-N or
%                   complex-valued analytic signal of length-2N)
%      tfd_params = cell of {win_length,win_type,[win_param],[doppler_or_time]}
%                   with:
%                     win_length  - length (odd value) of window,
%                     win_type    - name of window, e.g. 'hamm' for
%                                   Hamming window
%                     win_param   - if window has a parameter value,
%                                   e.g. 0.2 for Tukey window
%                     dopp_or_time - either:
%                       0 = to define window in the Doppler-domain
%                       1 = to define window in the time-domain
%      Ntime      = number of sample points in the time direction for
%                   the TFD, Ntime <= 2N
%
%  OUTPUT:
%     tf = Ntime x N TFD
%     R  = time-lag signal function
%     G1 = Doppler-independent kernel
%
%  EXAMPLE:
%     N=512; Ntime=128;
%     x=gen_LFM(N,0.1,0.3)+gen_LFM(N,0.4,0.1);
%
%     tf=dtfd_LI(x,{53,'hamm',0,1},Ntime); 
%     clf; vtfd(tf,x);


%
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
%
%--------------------------------------------------------------------------------
function [tf,R,g] = dtfd_LI(z,tfd_params,Ntime)
if(nargin<3) Ntime=[]; end
[z,N2,N,Nh]=check_analytic_signal(z);
tf=[]; R=[]; g2=[];


global_vars;
set_DBvars('all',0); set_DBvars('DBtfd',0);


%-------------------------------------------------------------------------
% 1. Get Doppler-window G_1[l]  (defined in the EITHER time domain 
%     or doppler domain) of length Q
%-------------------------------------------------------------------------
[G1,Q,time_doppler,Ntime]=get_Doppler_kernel(tfd_params,N,Ntime);
Qh=floor(Q/2);
Nh_time=floor(Ntime/2);
if(FAST_DTFD_DB) dispVars(Ntime,Q,Qh); end



%-------------------------------------------------------------------------
% 2. Multiple DF function K(l/2NT,k/2NT) by kernel G_1(l/2NT) 
%
%  extent: 0<nu<1/T;  g_1(t) is length Q.
%-------------------------------------------------------------------------
Z=fft(z);
R=zeros(Nh_time+1,N);
even_factor=rem(Q,2)-1;

l=0:Qh; 
for k=0:N-1
  i1=mod(k+l,N2); i2=mod(k-l,N2); 
  R(l+1,k+1) = Z(i1+1).*conj(Z(i2+1)).*G1(l+1);  
end

l=1:(Qh+even_factor); 
for k=0:N-1
  i1=mod(k+(N-l),N2); i2=mod(k-(N-l),N2); 
  R(Nh_time-l+1,k+1) = Z(i1+1).*conj(Z(i2+1)).*G1(Q-l+1);  
end

% Special case for l=Nh_time
k=0:N-1; l=Nh_time;
i1=mod(k+N,N2); i2=mod(k-N,N2); 
R(l+1,k+1) = Z(i1+1).*conj(Z(i2+1)).*G1(1);  


%-------------------------------------------------------------------------
% 3.  Expand R from:
%      0<nu<1/T
%    to:
%      |nu|<1/T 
%-------------------------------------------------------------------------
R=get_conj_dfR(R,Ntime,Nh_time,N);


%-------------------------------------------------------------------------
% 4. To time--frequency domain, using conjugate symmetry property.
%-------------------------------------------------------------------------
tf=ifft_conj_symm(R);
tf=tf./N2;  


tf=istfd_real(tf,z,mfilename);





function Rnew=get_conj_dfR(R,Ntime,Nh_time,N)
%--------------------------------------------------------------------------------
% Expand doppler--frequency function R[l,k]
%
%          from: 0<nu<1/T
%            to: |nu|<1/T
%
% R matrix has uniform (l/NT,k/2NT) discrete grid.
%--------------------------------------------------------------------------------
global_vars;
Rnew=zeros(Ntime,N);
if(FAST_DTFD_DBvars) dispVars(size(R),Nh_time,N); end
Rnew(1:Nh_time+1,:)=R(:,:);
even_factor=rem(Ntime,2)-1;


l=1:Nh_time+even_factor;
Rnew(Ntime-l+1,:)=conj( R(l+1,:) );
