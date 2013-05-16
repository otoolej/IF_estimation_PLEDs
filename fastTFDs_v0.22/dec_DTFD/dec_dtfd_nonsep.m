%-------------------------------------------------------------------------------
% Decimated TFD (time-frequency distribution)
% Kernel Type:  nonseparable
%
%
% USE: [tf,g] = dec_dtfd_nonsep(z,tfd_type,tfd_params,time_dec,freq_dec)
%
% INPUT: 
%      z = input signal (either real-valued signal of length-N or
%          complex-valued analytic signal of length-2N)
%
%      tfd_type = { 'wvd' | 'swvd' | 'pwvd' | 'sep' | 'cw' | 'mb' }
%            wvd  - Wigner-Ville distribution
%            swvd - Smoothed Wigner-Ville distribution
%                   (lag-independent kernel)
%            pwvd - Pseudo Wigner-Ville distribution
%                   (Doppler-independent kernel)
%            sep  - Separable-kernel distribution 
%                   (combintation of SWVD and PWVD)
%            cw   - Choi-Williams distribution
%            mb   - Modified-B distribution
% 
%      tfd_params = cell of kernel parameters:
%            wvd  - {}
%            swvd - {win_length,win_type,[win_param]}
%                   e.g. {11,'hamm'}
%            pwvd - {win_length,win_type,[win_param]}
%                   e.g. {200,'cosh',0.1
%            sep  - { {win1_length,win1_type,[win1_param]}, 
%                    {win2_length,win2_type,[win2_param]}
%                   where win1 is the doppler window and win2 is the 
%                   lag window, e.g. { {11,'hamm'}, {200,'cosh',0.1} }
%            cw   - {sigma_parameter}
%            mb   - {beta_parameter} in the range 1<beta<0
%
%      time_dec = integer value time decimation value
%                 where (2N/time_dec) must be an integer
%
%      freq_dec = integer value frequency decimation value
%                 where (N/freq_dec) must be an integer
%      
% OUTPUT: 
%      tf = LxJ TFD, where L=(2N/time_dec) and J=(N/freq_dec)
%      g  = kernel
%
% EXAMPLE:
%      N=10000; 
%      time_dec=160; freq_dec=80;
%      x=gen_LFM(N,0.1,0.3)+gen_LFM(N,0.4,0.1);
%
%      c=dec_dtfd_nonsep(x,'cw',{100},time_dec,freq_dec); 
%      vtfd(c,x);


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
function [tf,g] = dec_dtfd_nonsep(z,tfd_type,tfd_params,time_dec,freq_dec)
if(nargin<4) time_dec=1; end
if(nargin<5) freq_dec=1; end
[z,N2,N,Nh]=check_analytic_signal(z);
tf=[]; R=[]; g=[];


global_vars;
set_DBvars('all',0); set_DBvars('dbtfd',0);


%-------------------------------------------------------------------------------
% 1. Check decimation parameters
%-------------------------------------------------------------------------------
[L,time_dec]=checkParams_seq(time_dec,N2);
[J,freq_dec]=checkParams_seq(freq_dec,N);
Lh=ceil(L/2); Jh=ceil(J/2); J2=2*J;
if(FAST_DTFD_DBvars) dispVars(N2,N,time_dec,freq_dec,L,Lh,J,J2,Jh); end


% if product-kernel then generate window now:
G1=get_prod_kernel(tfd_type,tfd_params,N,N2);


%-------------------------------------------------------------------------
% 3. Get TL kernel Kdown[n,m]
%-------------------------------------------------------------------------
nt=0:time_dec:N2-1;
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
AFf_EvenM=zeros(Le,1);
AFf_OddM=zeros(Lo,1);


%-------------------------------------------------------------------------
% For m=even
%-------------------------------------------------------------------------
n=0:N-1; m1=1:Jh-1;
for m=0:Jh
  AFslice_EvenM=zeros(N,1);  
  
  %------------------------------------
  % Fold in the lag direction
  %------------------------------------
  for p=0:freq_dec-1
      g=gen_kern_lag_slice( 2*(p*J+m),tfd_type,tfd_params,N,G1);
      AFslice_EvenM = AFslice_EvenM + ...
          fft(getKslice_EvenM(z,n,p.*J+m,N2,N,Nh)) .* g;
  end

  %------------------------------------
  % Fold in the Doppler direction
  %------------------------------------
  AFf_EvenM=fold_sig_full(AFslice_EvenM,Le,time_dec_h);
  
  
  %------------------------------------
  % DFT the lag slice to the time--lag 
  % domain.
  %------------------------------------
  Racross(ni_even+1,m+1)=ifft_complex(AFf_EvenM);

end


%-------------------------------------------------------------------------
% For m=odd  (only if time_dec is odd value, otherwise no m=odd!)
%-------------------------------------------------------------------------
if(Lo>0)
  m2=0:Jh-1;
  if(FAST_DTFD_DBvars) dispVars((time_dec_h-1)/2); end
  for m=0:Jh
    AFslice_OddM=zeros(N,1);  
    
    %------------------------------------
    % Fold in the lag direction
    %------------------------------------
    for p=0:freq_dec-1
        g=gen_kern_lag_slice(2*(p*J+m)+1,tfd_type,tfd_params,N,G1);        
        AFslice_OddM = AFslice_OddM + ...
            fft(getKslice_OddM(z,n,p.*J+m,N2,N,Nh)) .* g;
    end

    %------------------------------------
    % Fold in the Doppler direction
    %------------------------------------
    if(time_dec>1)
      AFf_OddM=fold_sig_full_odd(AFslice_OddM,Lo,N,time_dec_h,...
                                 (time_dec_h-1)/2);  
    else
      AFf_OddM=AFslice_OddM;
    end
    

    %------------------------------------
    % DFT the lag slice to the time--lag 
    % domain.
    %------------------------------------
    Racross(ni_odd+1,m+1)=ifft_complex(AFf_OddM);
    
  end
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
tf=tf./N2;  


% final error check and return real TFD:
tf=istfd_real(tf,z,'dec_dtfd_nonsep');





function G1=get_prod_kernel(tfd_type,tfd_params,N,N2)
%---------------------------------------------------------------------
% Generate window for product-kernel
%---------------------------------------------------------------------

if( strcmp(tfd_type,'RID')==1 || strcmp(tfd_type,'prod')==1 || ...
    strcmp(tfd_type,'product')==1 )
    
    % oversample the window:
    L_dopp=N*N2;
    
    G1=get_Doppler_kernel(tfd_params,L_dopp,N2);
    G1=real(G1);
else
    G1=[];
end






function Kslice=getKslice(z,m,N,N2)
%-------------------------------------------------------------------------------
% Get TL function for one value of m, from n=0,...,N-1.
% (Kslice is Kdown[n,m_0])
%-------------------------------------------------------------------------------

Kslice=zeros(N,1);
n=0:N-1;
if(rem(m,2)==0)
  m=m/2;
  i1=mod(n+m,N2); i2=mod(n-m,N2);
  Kslice=z(i1+1).*conj(z(i2+1));
else
  m=(m-1)/2;
  i3=mod(n+m+1,N2); i2=mod(n-m,N2);
  Kslice = z(i3+1).*conj(z(i2+1)); 
end


function Kslice=getKslice_EvenM(z,n,m,N2,N,Nh)
%-------------------------------------------------------------------------------
% Get TL function for one value of m, from n=0,...,N-1.
% (Kslice is Kdown[n,m_0])
%-------------------------------------------------------------------------------
if(m<=Nh)
  i1=mod(n+m,N2); i2=mod(n-m,N2);
  Kslice=z(i1+1).*conj(z(i2+1));
else
  m=N-m;
  i1=mod(n+m,N2); i2=mod(n-m,N2);
  Kslice=conj(z(i1+1)).*z(i2+1);
end

  

function Kslice=getKslice_OddM(z,n,m,N2,N,Nh)
%-------------------------------------------------------------------------------
% Get TL function for one value of m, from n=0,...,N-1.
% (Kslice is Kdown[n,m_0])
%-------------------------------------------------------------------------------
if(m<Nh)
  i1=mod(n+m+1,N2); i2=mod(n-m,N2);
  Kslice=z(i1+1).*conj(z(i2+1));
else
  m=N-m-1;
  i1=mod(n+m+1,N2); i2=mod(n-m,N2);
  Kslice=conj(z(i1+1)).*z(i2+1);
end



