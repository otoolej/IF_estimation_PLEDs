%-------------------------------------------------------------------------------
% Decimated WVD (Wigner-Ville distribution).  See [1] for details.
%
%
% USE: wvd=dwvd_grid2(z,time_dec,freq_dec)
%
% INPUT:
%      z = N-point real-valued time-domain signal
%      time_dec  = either a constant, where U=(N/freq_dec) is an integer
%                  OR
%                  set freq_dec=[n1,n2,...,nU], where U<2N and 0<=ni<=2N 
%      freq_dec  = constant where J=(N/freq_dec) is an integer.
%
% OUTPUT:
%     wvd = U x J WVD, where
%                  - U is the length of the set time_dec (or if time_dec is a scalar
%                    then U=2N/time_dec)
%                  - J=(N/freq_dec)
%     K   = time-lag signal function
%
% EXAMPLE:
%      N=1024; freq_dec=4; 
%      time_dec=[10:5:1000,1000:1500,1700:10:2000];
%      x=gen_LFM(N,0.1,0.3)+gen_LFM(N,0.4,0.1);
%
%      wc=dwvd_grid2(x,time_dec,freq_dec); 
%      clf; vtfd(wc);

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
function [w,K] = dwvd_grid2(z,time_dec,freq_dec)
if(nargin<2) time_dec=1; end
if(nargin<3) freq_dec=1; end

[z,N2,N,Nh]=check_analytic_signal(z);
global_vars;
set_DBvars('all',0);  set_DBvars('dbtfd',0);
w=[]; K=[];



%-------------------------------------------------------------------------------
% 0. Check decimation parameters
%-------------------------------------------------------------------------------
n=checkParams_set(time_dec,N2);
[J,freq_dec]=checkParams_seq(freq_dec,N);
if(FAST_DTFD_DBvars) dispVars(time_dec,freq_dec,J); end


% First sort out even--odd time indices
ni_even=find( rem(n,2)==0 );
ni_odd=find( rem(n,2)==1 );


%-------------------------------------------------------------------------------
% 1. Get TL function for K(ntime,) folded in the m direction for
% frequency decimation.
%-------------------------------------------------------------------------------
K=get_Ktiaf(z,N2,N,Nh,n,ni_even-1,ni_odd-1,J,freq_dec); 


%-------------------------------------------------------------------------------
% 2. FFT to the TF domain.
%-------------------------------------------------------------------------------
w=dft_Racross(K,length(n),J,N,ni_even-1,ni_odd-1);
w=w./N2;  


w=istfd_real(w,z,'dwvd_grid2');






function K=get_Ktiaf(z,N2,N,Nh,n,ni_even,ni_odd,J,freq_dec)
%-------------------------------------------------------------------------
% 1. Get TL kernel K(nT/2,mT) for 0<t<2NT, 0<tau<(NT+1)
%    using K_across matrix (i.e. sequence points at (T/2,2T))
%
%    extent: |t|<NT and 0<tau<(NT+1)
%-------------------------------------------------------------------------
N_even=length(ni_even);
N_odd=length(ni_odd);
n_even=n(ni_even+1); n_odd=n(ni_odd+1);

K=zeros(N_even+N_odd,J);
Kslice=zeros(N,1);

N_odd_factor=rem(N,2); J_odd_factor=rem(J,2);
Jh=ceil(J/2);


% for n even values:
m=0:Nh; ji=0:Jh; mh=1:Nh-1; jih=1:Jh-1;
for in=1:N_even
  ne=n_even(in)/2;
  i1=mod(ne+m,N2); i2=mod(ne-m,N2);
  
  Kslice(m+1)=z(i1+1).*conj(z(i2+1)); 
  Kslice(N-mh+1)=conj( Kslice(mh+1));  
  
  % If need to fold for frequency decimation then do so here:
  if(J~=N)
    K(ni_even(in)+1,ji+1)=fold_sig_half(Kslice,J,Jh,freq_dec);
    K(ni_even(in)+1,J-jih+1)=conj( K(ni_even(in)+1,jih+1) );    
  else
    K(ni_even(in)+1,:)=Kslice;
  end
end


% for n odd values:
mh=0:Nh-1; jih=0:Jh-1;
for in=1:N_odd
  no=(n_odd(in)-1)/2;
  i3=mod(no+m+1,N2); i2=mod(no-m,N2);

  Kslice(m+1)=z(i3+1).*conj(z(i2+1)); 
  Kslice(N-mh)=conj(Kslice(mh+1));    

  % If need to fold for frequency decimation then do so here:
  if(J~=N)
    K(ni_odd(in)+1,ji+1)=fold_sig_half(Kslice,J,Jh,freq_dec);
    K(ni_odd(in)+1,J-1-jih+1)=conj(K(ni_odd(in)+1,jih+1));    
  else
    K(ni_odd(in)+1,:)=Kslice;
  end
end


