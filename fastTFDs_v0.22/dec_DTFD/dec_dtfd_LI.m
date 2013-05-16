%-------------------------------------------------------------------------------
% Decimated TFD (time-frequency distribution)
%
% Kernel Type: lag-independent kernel.  See [1] for details.
%
%
%  USE: [tf,R,G1] = dec_dtfd_LI(z,tfd_params,Ntime,time_dec,freq_dec)
%
%  INPUT:
%      z          = input signal (either real-valued signal of length-N or
%                   complex-valued analytic signal of length-2N)
%      tfd_params = cell of {win_length,win_type,[win_param],[dopp_or_time]}
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
%      time_dec   = constant where Ltime=(Ntime/time_dec) is an integer
%      freq_dec   = either a constant, where V=(N/freq_dec) is an integer
%                   OR
%                   freq_dec=[k1,k2,...,kV], where V<N and 0<=ki<=N 
%
%  OUTPUT:
%     tf = Ltime x V TFD, where
%                   - V is the length of set freq_dec (or V=N/freq_dec if
%                     freq_dec is a constant) 
%                   - Ltime=(Ntime/time_dec)
%     R  = time-lag signal function
%     G1 = lag-independent kernel
%
%  EXAMPLE:
%     t1=gen_LFM(171,0.1,0.4);
%     tf=dec_dtfd_LI(t1,{71,'hamm',[],1},512,4,[50:100,120:130]); 
%     vtfd(tf,t1);

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
function [tf,R,G1] = dec_dtfd_LI(z,tfd_params,Ntime,time_dec,freq_dec)
if(nargin<3) Ntime=[]; end
if(nargin<4) time_dec=1; end
if(nargin<5) freq_dec=1; end
[z,N2,N,Nh]=check_analytic_signal(z);
tf=[]; R=[]; g=[]; 


global_vars;
set_DBvars('all',0); set_DBvars('dbtfd',0);



%-------------------------------------------------------------------------
% 1. Get kernel g(nu,tau) = G_1(nu)  (defined in the EITHER time domain 
%     or doppler domain)
%
%  extent: 0<nu<1/T;  G_1(nu) is length Q.
%-------------------------------------------------------------------------
[G1,Q,time_doppler,Ntime]=get_Doppler_kernel(tfd_params,N,Ntime);
Qh=floor(Q/2);


%-------------------------------------------------------------------------------
% 2. Check decimation parameters
%-------------------------------------------------------------------------------
k=checkParams_set(freq_dec,N);
[L,time_dec]=checkParams_seq(time_dec,Ntime);
J=length(k); Lh=ceil(L/2);
if(FAST_DTFD_DBvars) dispVars(time_dec,freq_dec,L,Lh,J); end


Nh_time=floor(Ntime/2);

if(FAST_DTFD_DB) dispVars(Ntime,Q,Qh); end
g=G1;


%-------------------------------------------------------------------------
% 3. Multiple DF function K(l/2NT,k/2NT) by kernel G_1(l/2NT) 
%
%  extent: 0<nu<1/T;  g_1(t) is length Q.
%-------------------------------------------------------------------------
Z=fft(z);
Rslice=zeros(Ntime,1); R=zeros(L,J);
Q_even_factor=rem(Q,2)-1;
Nt_even_factor=rem(Ntime,2)-1;          

l=0:Qh;  l2=1:(Qh+Q_even_factor);  l3=1:(Nh_time+Nt_even_factor);
li=0:Lh; lbL=1:(Lh-1);
for ki=1:J
  
  %-------------------------------------------------------------------------
  % 4. From the frequency slice Kslice and multiply by the Doppler
  %    window G1[l]
  %-------------------------------------------------------------------------
  i1=mod(k(ki)+l,N2); i2=mod(k(ki)-l,N2); 
  Rslice(l+1) = Z(i1+1).*conj(Z(i2+1)).*G1(l+1);  

  i1=mod(k(ki)+(N-l2),N2); i2=mod(k(ki)-(N-l2),N2); 
  Rslice(Nh_time-l2+1) = Z(i1+1).*conj(Z(i2+1)).*G1(Q-l2+1); 
  
  i1=mod(k(ki)+N,N2); i2=mod(k(ki)-N,N2); 
  Rslice(Nh_time+1) = Z(i1+1).*conj(Z(i2+1)).*G1(1);  
  
  
  % get negative Doppler values from positive ones:
  Rslice(Ntime-l3+1)=conj( Rslice(l3+1) );

  %-------------------------------------------------------------------------
  % 5. If required, fold the frequency slice for time-domain decimation 
  %-------------------------------------------------------------------------
  if(L~=N2)
    R(li+1,ki)=fold_sig_half(Rslice,L,Lh,time_dec);
    R(li+1,ki)=R(li+1,ki)./time_dec;
    R(L-lbL+1,ki)=conj( R(lbL+1,ki) );          
  else
    R(:,ki)=Rslice;
  end
  
end


tf=ifft_conj_symm(R);
tf=tf./Ntime;  


tf=istfd_real(tf,z,mfilename);





