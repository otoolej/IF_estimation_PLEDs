%-------------------------------------------------------------------------------
% Decimated TFD (time-frequency distribution)
%
% Kernel Type: separable-kernel.  See [1] for details.
%
% 
% USE: [tf,g] = dec_dtfd_sep(z,doppler_win_params,lag_win_params,...
%                                  Ntime,Nfreq,time_dec,freq_dec)
%
% INPUT: 
%      z = input signal (either real-valued signal of length-N or
%          complex-valued analytic signal of length-2N)
%
%      doppler_win_params = {win_length,win_type,[win_param],win_DFT} 
%                   with:
%                       win_length - length (odd value) of window,
%                       win_type   - name of window, e.g. 'hamm'
%                                    for Hamming window
%                       win_param  - if window has a parameter
%                                    value, e.g. 0.2 for Tukey window
%                       win_DFT    - to use the DFT of window,
%                                    0=no, 1=yes (default=0)
%                   e.g. {11,'hamm',0,1} or {200,'cosh',0.1} 
%
%      lag_win_params = {win_length,win_type,[win_param],win_DFT} 
%                   e.g. {111,'hamm'} or {127,'tukey',0.2} 
%
%      Ntime = number of sample points in the time direction for 
%              the TFD, Ntime <= 2N
%      Nfreq = number of sample points in the frequency direction for
%              the TFD, Nfreq <= N
%      time_dec = integer value time decimation value
%                 where (Ntime/time_dec) must be an integer
%      freq_dec = integer value frequency decimation value
%                 where (Nfreq/freq_dec) must be an integer
%      
%
% OUTPUT: 
%      tf = Ltime x Jfreq TFD, where
%                - Ltime=(Ntime/time_dec)
%                - Jfreq=(Nfreq/freq_dec)
%      g  = Doppler-lag kernel 
%
%
% EXAMPLE:
%      N=10000; Ntime=512; Nfreq=256;
%      time_dec=4; freq_dec=2;
%      x=gen_LFM(N,0.1,0.3)+gen_LFM(N,0.4,0.1);
%
%      c=dec_dtfd_sep(x,{51,'hamm',0,1},{271,'hann'},Ntime,Nfreq,time_dec,freq_dec); 
%      vtfd(c,x);
    

%
%  Copyright (c) 2010 John M. O' Toole, The University of Queensland
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
function [tf,g]=dec_dtfd_sep(z,doppler_win_params,lag_win_params, ...
                             Ntime,Nfreq,time_dec,freq_dec)
if(nargin<4) Ntime=[]; end
if(nargin<5) Nfreq=[]; end
if(nargin<6) time_dec=1; end
if(nargin<7) freq_dec=1; end
[z,N2,N,Nh]=check_analytic_signal(z);
tf=[]; R=[]; g=[];


global_vars;
DB_FOR_THIS_FILE=0;
set_DBvars('all',DB_FOR_THIS_FILE); 
set_DBvars('dbtest',DB_FOR_THIS_FILE);



%-------------------------------------------------------------------------
% 0.1 Get doppler window G_1[l]  (defined in the EITHER time domain 
%     or doppler domain)
%-------------------------------------------------------------------------
[G1,Q,time_doppler,Ntime]=get_Doppler_kernel(doppler_win_params,N, ...
                                                          Ntime);
% If Ntime is odd then adjust and check out G1 again: (should
% really do this check in 'get_Doppler_kernel.m' function)
if(rem(Ntime,2)==1)
  Ntime=Ntime+1;
  [G1,Q,time_doppler,Ntime]=get_Doppler_kernel(doppler_win_params,N, ...
                                                          Ntime);
end

Qh=floor(Q/2);
Nh_time=floor(Ntime/2);
Nhtime_even_factor=rem(Ntime,2)-1;
Q_even_factor=rem(Q,2)-1;
if(FAST_DTFD_DBvars)
  dispVars(Ntime,Nh_time,Nhtime_even_factor,Q,Qh,Q_even_factor); 
end
g=G1;
G1_pad=padWin(G1,Nh_time);


if(Nh_time~=N)
  % Modulate signal for ood m values 
  G1mod=zeros(Q,1);
  lp=0:Qh; ln=1:(Qh+Q_even_factor);
  G1mod(lp+1)=G1(lp+1).*exp(-j*pi/N.*lp(:)).*exp(j*pi/Nh_time.*lp(:));
  G1mod(Q-ln+1)=G1(Q-ln+1).*exp(-j*pi/N.*(N-ln(:))).* ...
      exp(j*pi/Nh_time.*(Nh_time-ln(:)));
else
  G1mod=G1;
end


%-------------------------------------------------------------------------
% 0.2 Get lag window g_2[m]
%-------------------------------------------------------------------------
[g2,P,Ph,Nfreq]=get_lag_kernel(lag_win_params,N2,N,Nfreq);
g2=g2.';
N2_freq=2*Nfreq;
Pq=ceil(Ph/2);
Ph_even_factor=rem(Ph,2)-1;
g2_pad_t=padWin(g2,N2_freq);
g2_pad=g2_padWin(g2,N2_freq,P,Ph);
Nh_freq=floor(Nfreq/2);

if(FAST_DTFD_DBvars) dispVars(Nfreq,P,Ph,Pq,Ph_even_factor); end



%-------------------------------------------------------------------------------
% 0.3 Check decimation parameters
%-------------------------------------------------------------------------------
[L,time_dec]=checkParams_seq(time_dec,Ntime);
[J,freq_dec]=checkParams_seq(freq_dec,Nfreq);
Lh=ceil(L/2); Jh=ceil(J/2); J2=2*J;
if(FAST_DTFD_DBvars) 
  dispVars(Ntime,Nfreq,time_dec,freq_dec,L,J,Jh); 
end




%-------------------------------------------------------------------------
% 0.4 Set up some variables 
%-------------------------------------------------------------------------
nt=0:time_dec:Ntime-1;
ni_even=find(rem(nt,2)==0)-1;
ni_odd=find(rem(nt,2)==1)-1;
Le=length(ni_even); Lo=length(ni_odd);

if(rem(time_dec,2)==0)
  time_dec_h=time_dec/2;  
else
  time_dec_h=time_dec;  
end
if(FAST_DTFD_DBvars) dispVars(Le,Lo,time_dec_h); end

Racross=zeros(Le+Lo,J);
rf_EvenM=zeros(Le,1);
rf_OddM=zeros(Lo,1);


%-------------------------------------------------------------------------
% For m=even
%-------------------------------------------------------------------------
n=0:N-1;
lp=0:Qh; ln=1:(Qh+Q_even_factor);



for m=0:(Pq-1-Ph_even_factor)
  Kslice_EvenM=zeros(N,1);  
  rslice_EvenM=zeros(Nh_time+1+Nhtime_even_factor,1);    
  
  %-------------------------------------------------------------------------
  % a) Form lag-slice of windowed time-lag function and then fold in
  %    the lag direction.
  %-------------------------------------------------------------------------
  for p=0:freq_dec-1
    Kslice_EvenM=Kslice_EvenM+ ... 
        getRtimeslice_EvenM(z,n,p.*J+m,N2,Nfreq,Nh_freq,g2_pad);
  end


  %-------------------------------------------------------------------------
  % b) DFT to the Doppler--lag domain..
  %-------------------------------------------------------------------------
  Kslice_EvenM=fft(Kslice_EvenM);

  %-------------------------------------------------------------------------
  % c) Multiply by Doppler window G1
  %-------------------------------------------------------------------------
  rslice_EvenM(lp+1)=Kslice_EvenM(lp+1).*G1(lp+1);
  rslice_EvenM(Nh_time-ln+1)=Kslice_EvenM(N-ln+1).* G1(Q-ln+1);  
  
  
  %-------------------------------------------------------------------------
  % d) Fold in the Doppler direction
  %-------------------------------------------------------------------------
  rf_EvenM=fold_sig_full(rslice_EvenM,Le,time_dec_h);
  
  
  %-------------------------------------------------------------------------
  % e) DFT the lag slice to the time--lag domain
  %-------------------------------------------------------------------------
  Racross(ni_even+1,m+1)=ifft_complex(rf_EvenM);
end

%-------------------------------------------------------------------------
% For m=odd  (only if time_dec is odd value, otherwise no m=odd!)
%-------------------------------------------------------------------------
if(Lo>0)
  if(FAST_DTFD_DBvars) dispVars((time_dec_h-1)/2); end
  for m=0:Pq-1
    Kslice_OddM=zeros(N,1);  
    rslice_OddM=zeros(Nh_time,1);    

    %-------------------------------------------------------------------------
    % a) Form windowed time-lag function and then fold in the lag direction.
    %-------------------------------------------------------------------------
    for p=0:freq_dec-1
      Kslice_OddM=Kslice_OddM+ ...
          getRtimeslice_OddM(z,n,p.*J+m,N2,Nfreq,Nh_freq,g2_pad);      
    end

    %-------------------------------------------------------------------------
    % b) DFT to the Doppler--lag domain.
    %-------------------------------------------------------------------------
    Kslice_OddM=fft(Kslice_OddM);

    %-------------------------------------------------------------------------
    % c) and multiply by Doppler window G1
    %-------------------------------------------------------------------------
    rslice_OddM(lp+1)=Kslice_OddM(lp+1).*G1mod(lp+1);
    rslice_OddM(Nh_time-ln+1)=Kslice_OddM(N-ln+1).* G1mod(Q-ln+1);  
    
  
    %-------------------------------------------------------------------------
    % d) Fold in the Doppler direction
    %-------------------------------------------------------------------------
    if(time_dec>1)
      rf_OddM=fold_sig_full_odd(rslice_OddM,Lo,N,time_dec_h, ...
                                (time_dec_h-1)/2);  
    else
      rf_OddM=rslice_OddM;
    end
  
  
    %-------------------------------------------------------------------------
    % e) DFT the lag slice to the time--lag domain.
    %-------------------------------------------------------------------------      
    Racross(ni_odd+1,m+1)=ifft_complex(rf_OddM);
  end
end


%-------------------------------------------------------------------------      
% JUST TO CHECK:
%-------------------------------------------------------------------------      
for m=(Pq-Ph_even_factor):Jh
  Racross(ni_even+1,m+1)=0;
end
for m=Pq:(Jh-1)
  Racross(ni_odd+1,m+1)=0;
end


%-------------------------------------------------------------------------
% Then recover negative lag values from positive ones
%-------------------------------------------------------------------------
m1=1:Jh-1; m2=0:Jh-1;
Racross(ni_even+1,J-m1+1)=conj( Racross(ni_even+1,m1+1) );
if(~isempty(ni_odd))
  Racross(ni_odd+1,J-m2)=conj( Racross(ni_odd+1,m2+1) );
end



tf=dft_Racross(Racross,L,J,N,ni_even,ni_odd);

tf=tf./Ntime;  
tf=istfd_real(tf,z,'dec_dtfd_sep');






function Rslice=getRtimeslice_EvenM(z,n,m,N2,Nfreq,Nh_freq,g2)
%-------------------------------------------------------------------------------
% Get TL function for one value of m, from n=0,...,N-1.
% (Kslice is Kdown[n,m_0])
%-------------------------------------------------------------------------------

if(m<=Nh_freq)
  i1=mod(n+m,N2); i2=mod(n-m,N2);
  Rslice=z(i1+1).*conj(z(i2+1)).*g2(2.*m+1);
else
  m1=Nfreq-m;
  i1=mod(n+m1,N2); i2=mod(n-m1,N2); ig=mod(2.*m1,Nfreq);
  Rslice=conj(z(i1+1)).*z(i2+1).*g2(ig+1);
end


function Rslice=getRtimeslice_OddM(z,n,m,N2,Nfreq,Nh_freq,g2)
%-------------------------------------------------------------------------------
% Get TL function for one value of m, from n=0,...,N-1.
% (Kslice is Kdown[n,m_0])
%-------------------------------------------------------------------------------

if(m<=Nh_freq-1)
  i1=mod(n+m+1,N2); i2=mod(n-m,N2); 
  Rslice=z(i1+1).*conj(z(i2+1)).*g2(2.*m+2);
else
  m1=Nfreq-m-1;
  i1=mod(n+m1+1,N2); i2=mod(n-m1,N2); ig=mod(2.*m1+1,2*Nfreq);
  Rslice=conj(z(i1+1)).*z(i2+1).*g2(ig+1);
end




function g2_pad=g2_padWin(g2,N2_freq,P,Ph)
%-------------------------------------------------------------------------------
% Pad window g2[m] from length P to 2*Nfreq
%-------------------------------------------------------------------------------
P_even_factor=rem(P,2)-1;
g2_pad=zeros(N2_freq,1);

m=0:Ph;
g2_pad(m+1)=g2(m+1);
m=1:(Ph+P_even_factor);
g2_pad(N2_freq-m+1)=g2(P-m+1);
