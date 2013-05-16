%-----------------------------------------------------------------------------------------
%  FFT time--lag function to time--frequency domain. 
%
%  USE:  tfd=dft_Rdown(R,Ntime,Nfreq)
%
%  INPUT:
%        R      - R_down matrix R(T,T)
%        Ntime  - number of sample points in the time direction
%        Nfreq  - number of sample points in the frequency direction
%
%  OUTPUT:
%       TFD array
%

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
%-----------------------------------------------------------------------------------------
function tfd=dft_Rdown(R,Ntime,Nfreq)
tfd=zeros(Ntime,Nfreq);
N2_freq=2*Nfreq;
n_even=0:2:Ntime-1; n_odd=1:2:Ntime-1;
m_even=0:2:N2_freq-1; m_odd=1:2:N2_freq-1;

DBprofile=0;

%---------------------------------------------------------------------
% 1. Even n values.
%---------------------------------------------------------------------
if(DBprofile) dispVars('DFT for n even.'); end
tfd(n_even+1,:)=fft_conj_symm(R(:,m_even+1).').';

%---------------------------------------------------------------------
% 2. Odd n values.  (time--lag function NOT conjugate symmetrical.)
%---------------------------------------------------------------------

%-------------------------------------------------
% 2a. for k=1,...,N-1  (i.e. not k=0)
%-------------------------------------------------
if(DBprofile) dispVars('Setup for n odd.'); end
n_odd=1:2:Ntime-1; k=1:Nfreq-1;
A=cos( (pi/Nfreq).*k ); B=sin( (pi/Nfreq).*k );
C=( (A.^2+B.^2)./B);


Rfold=get_imagDFT_part(R(:,m_odd+1));
if(DBprofile) dispVars('DFT for n odd.'); end
tfd_part=fft_conj_symm(Rfold.').';

n=0;
for n_odd=1:2:Ntime-1
  tfd(n_odd+1,k+1)=tfd_part(n+1,k+1).*C;
  n=n+1;
end

%-------------------------------------------------
% 2b. for k=0;
%-------------------------------------------------
n=0;
for n_odd=1:2:Ntime-1
  tfd(n_odd+1,1)=sum( R(n+1,m_odd+1) );
  n=n+1;
end



%---------------------------------------------------------------------
% 
%---------------------------------------------------------------------
function Rfold=get_imagDFT_part(R)
[N,Nfreq]=size(R);  Nh=ceil(N/2); Nh_freq=ceil(Nfreq/2);
n=0:N-1;
Rfold=zeros(N,Nfreq);


Rfold(n+1,1)=imag(R(n+1,1));

k=1:Nh_freq; ke=(Nh_freq+1):Nfreq-1;
for n=0:N-1
  Rfold(n+1,k+1)=(R(n+1,k+1) - conj(R(n+1,Nfreq-k+1)))./(2*j);
  Rfold(n+1,ke+1)=conj(Rfold(n+1,Nfreq-ke+1));
end




