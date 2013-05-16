%--------------------------------------------------------------------------------
% Get lag kernel G_1(l/(QT))
%
%  USE:  [G1,G1_length,time_doppler,Ntime]=get_Doppler_kernel(tfd_params,N,Ntime)
%
%  INPUT:
%       tfd_params - cell of {win_length,win_type,[win_param],[Doppler_or_time]}
%                     e.g. {11,'hamm',0,1}
%                     the parameter Doppler_or_time is either:
%                       0 = to define window in the time-domain
%                       1 = to define window in the Doppler-domain
%       N          - length of signal
%       Ntime      - number of sample points in the time direction
%
%  OUTPUT:
%       G1           - Doppler function G_1(1/(QT))
%       G1_length    - length of G1
%       time_doppler - string either of the value 'time' or 'doppler' 
%       Ntime        - returns Ntime as it may be adjusted from input
%       
% EXAMPLE:
%       G1=get_Doppler_kernel({111,'cosh',0.1,1},128,512);
%       plot(real(G1));
%
% 
% IF window is defined in doppler domain:
%     window function G_1(l/NT) is length Q
% OR if window is defined in time domain:
%     window function g_1(nT) is length-Q, and then zero-padded 
%     to the doppler domain where G_1(l/NT) is length N
% ENDIF 
%
% Also, sets an appropriate value for Ntime if not already set.
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
%
%--------------------------------------------------------------------------------
function [G1,G1_length,time_doppler,Ntime]=get_Doppler_kernel(tfd_params,N,Ntime)
if(nargin<3) error('Need G1 parameters and signal length'); end

% If calling from outside then maybe do not need to set this here:
global_vars; 
set_DBvars('all',0); set_DBvars('dbwarn',0); 
set_DBvars('dbplot',0); set_DBvars('dbvars',0);



l=length(tfd_params);
if( l<2 ) 
  error('Need at least two window parameters'); 
end
Q=tfd_params{1};
% not necessary but good to force this:
if(Q/2==fix(Q/2)) 
  Q=Q-1;
  if(FAST_DTFD_DBwarn)
    dispDB(['Forcing Q to be odd. Q is now Q=',num2str(Q)]); 
  end
end

win_type=tfd_params{2}; 
win_param=0; win_param2=0;
if( l>=3 ) win_param=tfd_params{3}; end
if( l>=4 ) win_param2=tfd_params{4}; end
Qh=ceil(Q/2);   Qh_floor=floor(Q/2);


%-------------------------------------------------------------------------
% If defining window in the Doppler domain. 
%-------------------------------------------------------------------------
if( win_param2==1 )
  if(isempty(Ntime)) 
    if(Q>N) Ntime=2*N;  else   Ntime=2*Q; end
  end
  Nh_time=floor(Ntime/2);
    
  if(Q>Nh_time)
    warning('Chopping down length of G1 to Nh_time');
    Q=make_odd(Nh_time-1);
  end
  g1=getWin(Q,win_type,win_param);

  time_doppler='doppler';
  G1=shiftWin(g1);
  G1pad=padWin(G1,Nh_time);

  
  %---------------------------------------------------------------------
  % If Nh_time>Q and Q even then window will be different
  % (see padWin.m for details).
  %---------------------------------------------------------------------
  if(Nh_time>Q && rem(Q,2)==0)
    Q=Q+1;
  end
  G1=zeros(Q,1); even_factor=rem(Q,2)-1;
  Qh_floor=floor(Q/2);
  
  l=0:Qh_floor;
  G1(l+1)=G1pad(l+1);
  l=1:Qh_floor+even_factor;
  G1(Q-l+1)=G1pad(Nh_time-l+1);


  if(FAST_DTFD_DBvars) dispVars(Ntime,Nh_time,Q); end
  
  if( isreal_fn(ifft(G1pad))==0 )
    disp('-- WARNING: g1 not real --');
  end
  G1_length=Q;
  
  
%-------------------------------------------------------------------------  
% or if defining window in time domain
%-------------------------------------------------------------------------
else
  time_doppler='time';
  if(Q>N)
    warning('Chopping down length of G1 to N');
    Q=N;
  end

  g1=getWin(Q,win_type,win_param);
  g1=shiftWin(g1);
  g1=padWin(g1,N);
  G1_length=N;


  G1=fft(g1); 
  G1=G1./G1(1);

  
  if( isreal_fn(g1)==0 )
    error('Window function g1(t) must be real valued.');
  end
  
  % If window is defined in the time domain then no time-decimation is allowed,
  % that is, Ntime should equal 2N.
  if(Ntime~=2*N | isempty(Ntime) ) Ntime=2*N; end
end

if(FAST_DTFD_DBwarn) 
  dispDB(['||- Doppler window defined in the ' time_doppler ' domain. -||']); 
end


if(FAST_DTFD_DBplot)
  figure(555); clf;
  plotc(ifft(G1(:)),1,'Time-window g1[n]'); 
end
