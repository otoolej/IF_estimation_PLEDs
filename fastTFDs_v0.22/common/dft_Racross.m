%-----------------------------------------------------------------------------------------
%  FFT time--lag function to time--frequency domain. 
%
%  USE:  tfd=dft_Racross(R,N2,Nfreq,N,n_even,n_odd)
%
%  INPUT:
%        R      - R_acroos matrix R(T/2,2T)
%        N2     - 2*N
%        Nfreq  - number of sample points in the frequency direction
%        N      - length of signal
%        n_even - set of even time-sample points
%        n_odd  - set of odd time-sample points
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
function tfd=dft_Racross(R,N2,Nfreq,N,n_even,n_odd)
if(nargin<5) n_even=0:2:N2-1; end
if(nargin<6) n_odd=1:2:N2-1; end
tfd=zeros(N2,Nfreq);


%---------------------------------------------------------------------
% 1. Even n values.
%---------------------------------------------------------------------
tfd(n_even+1,:)=fft_conj_symm(R(n_even+1,:).').';


%---------------------------------------------------------------------
% 2. Odd n values.  (time--lag function NOT conjugate symmetrical.)
%---------------------------------------------------------------------

%-------------------------------------------------
% 2a. for k=1,...,N-1  (i.e. not k=0)
%-------------------------------------------------
k=1:Nfreq-1;
A=cos( (pi/Nfreq).*k ); B=sin( (pi/Nfreq).*k );
C=((A.^2)./B+B);

Rfold=get_imagDFT_part(R(n_odd+1,:));
tfd_part=fft_conj_symm(Rfold.').';

for ni=1:length(n_odd)
  tfd(n_odd(ni)+1,k+1)=tfd_part(ni,k+1).*C;
end

%-------------------------------------------------
% 2b. for k=0;
%-------------------------------------------------
for ni=1:length(n_odd)
  tfd(n_odd(ni)+1,1)=sum( R(n_odd(ni)+1,:) );
end






function Rfold=get_imagDFT_part(R)
[Nt,Nfreq]=size(R); Nh_freq=ceil(Nfreq/2);
n=0:Nt-1;
Rfold=zeros(Nt,Nfreq);


Rfold(n+1,1)=imag(R(n+1,1));

k=1:Nh_freq; ke=(Nh_freq+1):Nfreq-1;
for n=0:Nt-1
  Rfold(n+1,k+1)=(R(n+1,k+1) - conj(R(n+1,Nfreq-k+1)))./(2*j);
  Rfold(n+1,ke+1)=conj(Rfold(n+1,Nfreq-ke+1));
end





