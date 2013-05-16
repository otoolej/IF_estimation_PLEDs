%--------------------------------------------------------------------------------
% 
%  Get lag kernel g_2(2mT)
%
%  USE: [g2,P,Ph_floor,Nfreq]=get_lag_kernel(win_params,N2,N,Nfreq)
%
%  INPUT:
%       win_params - cell of {win_length,win_type,[win_param],[lag_or_freq]}
%                     e.g. {11,'hamm',0,1}
%                     the parameter lag_or_freq is either:
%                       0 = to define window in the lag-domain
%                       1 = to define window in the frequency-domain
%       N2         - 2*N
%       N          - lenght of signal
%       Nfreq      - number of sample points in the frequency direction
%
%  OUTPUT:
%       g2       - lag kernel g_2(2mT)
%       P        - length of g2 (as may have changed)
%       Ph_floor - floor(P/2)
%       Nfreq    - Nfreq (as may have changed)
%
%  EXAMPLE:
%       g2=get_lag_kernel({61,'hann'},256,128,512);
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
%--------------------------------------------------------------------------------
function [g2,P,Ph_floor,Nfreq,g2_pad]=get_lag_kernel(win_params,N2,N,Nfreq)
if(nargin<2) error('Need Nfreq value too'); end

% If calling from outside then maybe do not need to set this here:
global_vars; 
set_DBvars('all',0); set_DBvars('dbwarn',0); 
set_DBvars('dbplot',0); set_DBvars('dbvars',0);


l=length(win_params);
if( l<2 ) 
  error('Need at least two window parameters'); 
end
P=win_params{1};
if(P>N2)
  warning('Chopping down length of g2 to 2N');
  P=N2;
end
% FORCE P to be odd.
if(P/2==fix(P/2)) 
  P=P-1;
  if(FAST_DTFD_DBwarn)
    dispDB(['Forcing P to be odd. P is now P=',num2str(P)]); 
  end
end

  

win_type=win_params{2}; 
win_extra_param=0;
win_dft_param=0;
if( l>=3 ) win_extra_param=win_params{3}; end
if( l>3 )  win_dft_param=win_params{4}; end

%---------------------------------------------------------------------
% and/or if reversing the window (in time)
%---------------------------------------------------------------------
reverse=0;
if(length(win_params)>5)
  reverse=win_params{6};
end



Ph_floor=floor(P/2);
Ph_ceil=ceil(P/2);
if(FAST_DTFD_DB) dispVars(P,win_type,win_extra_param,win_dft_param); end


g2=getWin(P,win_type,win_extra_param,win_dft_param);
g2=g2(:);

%---------------------------------------------------------------------
% Shift to appropriate indicies 
%---------------------------------------------------------------------
g2=shiftWin(g2);


%---------------------------------------------------------------------
% Work out an appropriate size for Nfreq
%---------------------------------------------------------------------
Ph_ceil_plus=ceil((P+1)/2);
if(isempty(Nfreq))
  if(P~=N2)
    Nfreq=Ph_ceil_plus;
  else
    Nfreq=N;
  end
end  

if( P~=N2 && Nfreq<Ph_ceil_plus  ) 
  Nfreq=Ph_ceil_plus;
  if(FAST_DTFD_DBwarn) dispDB('Nfreq too small. Increasing.'); end
end
if( P==N2 && Nfreq<N )
  Nfreq=N;
  if(FAST_DTFD_DBwarn) dispDB('Nfreq too small. Increasing.'); end
end



N2freq=2*Nfreq;
g2_pad=padWin(g2,N2freq);


%---------------------------------------------------------------------
% If 2*Nfreq>P and P even then g2 and g2_pad windows will be different
% (see padWin.m for details).
%---------------------------------------------------------------------
if(N2freq>P && rem(P,2)==0 && P~=N2)
  P=P+1;
  if(FAST_DTFD_DBwarn) 
    dispDB(['Increasing (zero-padding) size of P to: P=',num2str(P)]); 
  end
  
  % need to return same window (g2=g2_pad in nonzero region):
  g2=zeros(P,1); even_factor=rem(P,2)-1;
  Ph_floor=floor(P/2);

  l=0:Ph_floor;
  g2(l+1)=g2_pad(l+1);
  l=1:Ph_floor+even_factor;
  g2(P-l+1)=g2_pad(N2freq-l+1);
end


%---------------------------------------------------------------------
% if normalizing window:
%---------------------------------------------------------------------
if(length(win_params)>4)
  if(strcmpi(win_params{5},'y')==1)
    g2=g2./sum(g2);
  end
end

%---------------------------------------------------------------------
% if (time) reversing window
%---------------------------------------------------------------------
if(reverse)
  g2=shiftWin(g2);
  g2_pad=shiftWin(g2_pad);
end


if( isreal_fn(fft(g2_pad))==0 )
  error('Window function G2(f) must be real valued.');
end



if(FAST_DTFD_DBplot)
  figure(333);
  subplot(2,1,1); plot(g2(1:Ph_floor+1));
  title('lag-window g2[m]');
  subplot(2,1,2); plot(g2_pad(1:Ph_floor+1));
  title('lag-window g2[m]');
  dispEE(g2(1:Ph_floor+1),g2_pad(1:Ph_floor+1))
end



if(FAST_DTFD_DBvars) dispVars(P,Ph_floor,Nfreq,N2freq); end
