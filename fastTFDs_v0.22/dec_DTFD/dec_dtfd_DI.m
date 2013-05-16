%-------------------------------------------------------------------------------
% Decimated TFD (time-frequency distribution)
%
% Kernel Type: Doppler-independent kernel.  See [1] for details.
%
%
%  USE: [tf,R,g2]=dec_dtfd_DI(z,tfd_params,Nfreq,time_dec,freq_dec)
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
%      time_dec   = either a constant, where U=(2N/time_dec) is an integer
%                   OR
%                   time_dec=[n1,n2,...,nU], where U<2N  and 0<=ni<=2N 
%      freq_dec   = constant where Jfreq=(N/freq_dec) is an integer
%
%  OUTPUT:
%     tf = U x Jfreq TFD, where
%                   - U is the length of set time_dec (or U=N2/time_dec if
%                     time_dec is a constant) 
%                   - Jfreq=(N/freq_dec)
%     R  = time-lag signal function
%     g2 = Doppler-independent kernel
%
%  EXAMPLE:
%     t1=gen_LFM(171,0.1,0.4);
%     tf=dec_dtfd_DI(t1,{51,'hamm'},258,[10:60,141:200,210,211],3); 
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
function [tf,R,g2] = dec_dtfd_DI(z,tfd_params,Nfreq,time_dec,freq_dec)
if(nargin<3) Nfreq=[]; end
if(nargin<4) time_dec=1; end
if(nargin<5) freq_dec=1; end

[z,N2,N,Nh]=check_analytic_signal(z);

tf=[]; R=[]; g2=[];


global_vars;
DB_FOR_FILE=0;
set_DBvars('all',DB_FOR_FILE); set_DBvars('dbtfd',DB_FOR_FILE);

%-------------------------------------------------------------------------
% 1. Get lag-window g_2[m]
%-------------------------------------------------------------------------
[g2,P,Ph,Nfreq]=get_lag_kernel(tfd_params,N2,N,Nfreq);
g2=g2(:);
Nh_freq=floor(Nfreq/2);
Pq=ceil(Ph/2);
if(FAST_DTFD_DBvars) dispVars(P,Ph,Pq,Nfreq,Nh_freq); end


%-------------------------------------------------------------------------------
% 2. Check decimation parameters
%-------------------------------------------------------------------------------
n=checkParams_set(time_dec,N2);
[J,freq_dec]=checkParams_seq(freq_dec,Nfreq);
L=length(n); Jh=ceil(J/2);
if(FAST_DTFD_DBvars) dispVars(time_dec,freq_dec,L,J,Jh); end


%-------------------------------------------------------------------------------
% 3. Get smoothed TL function for K(ntime,) folded in the m direction for
%    frequency decimation.
%-------------------------------------------------------------------------------
% First sort out even--odd time indices
ni_even=find( rem(n,2)==0 )-1;
ni_odd=find( rem(n,2)==1 )-1;
n_even=n(ni_even+1); n_odd=n(ni_odd+1);
N_even=length(ni_even);
N_odd=length(ni_odd);

R=zeros(N_even+N_odd,J);
Rslice=zeros(Nfreq,1);

Ph_even_factor=rem(Ph,2)-1;

%------------------------------------
% for n even values:
%------------------------------------
m=0:(Pq-1-Ph_even_factor);
ji=0:Jh; mh=1:Nh_freq-1; jih=1:Jh-1;
for in=1:N_even
  ne=n_even(in)/2;
  i1=mod(ne+m,N2); i2=mod(ne-m,N2);
  
  Rslice(m+1)=z(i1+1).*conj(z(i2+1)).*g2(2.*m+1); 
  Rslice(Nfreq-mh+1)=conj( Rslice(mh+1));  
  
  % If need to fold for frequency decimation then do so here:
  if(J~=N)
    R(ni_even(in)+1,ji+1)=fold_sig_half(Rslice,J,Jh,freq_dec);
    R(ni_even(in)+1,J-jih+1)=conj( R(ni_even(in)+1,jih+1) );    
  else
    R(ni_even(in)+1,:)=Rslice;
  end
end


%------------------------------------
% for n odd values:
%------------------------------------
m=0:(Pq-1);
mh=0:Nh_freq-1; jih=0:Jh-1;
for in=1:N_odd
  no=(n_odd(in)-1)/2;
  i3=mod(no+m+1,N2); i2=mod(no-m,N2);

  Rslice(m+1)=z(i3+1).*conj(z(i2+1)).*g2(2.*m+2); 
  Rslice(Nfreq-mh)=conj(Rslice(mh+1));    

  % If need to fold for frequency decimation then do so here:
  if(J~=N)
    R(ni_odd(in)+1,ji+1)=fold_sig_half(Rslice,J,Jh,freq_dec);
    R(ni_odd(in)+1,J-1-jih+1)=conj(R(ni_odd(in)+1,jih+1));    
  else
    R(ni_odd(in)+1,:)=Rslice;
  end
end


%-------------------------------------------------------------------------------
% 3. FFT to the TF domain.
%-------------------------------------------------------------------------------
tf=dft_Racross(R,length(n),J,N,ni_even,ni_odd);
tf=tf./N2;  


tf=istfd_real(tf,z,'dec_dtfd_DI');

