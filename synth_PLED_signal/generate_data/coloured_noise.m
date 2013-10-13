%-------------------------------------------------------------------------------
% coloured_noise: generate either white or coloured Guassian noise (see [1])
%
% Syntax: n=coloured_noise(TYPE,N,Fs)
%
% Inputs: 
%     TYPE - type of noise, either WGN or CGN (white or coloured Gaussian noise)
%     N    - length of signal
%     Fs   - sample frequency
%
% Outputs: 
%     n    - noise
%
% Example:
%     N=1024; Fs=50;
%     ns=coloured_noise('CGN',N,Fs);
%    
%     figure(1); clf;
%     subplot(2,1,1); n=(0:N-1)./Fs
%     plot(n,ns); xlabel('time (seconds)');
%     subplot(2,1,2); k=(0:N-1)./(N/Fs);
%     plot( k,abs(fft(ns)).^2 ); xlim(xlim./2); xlabel('frequency (Hz)');
%
%
% [1] Kasdin, N. Jeremy. "Discrete simulation of colored noise and stochastic processes and
% 1/f^α power law noise generation." Proc. IEEE, 83(5), 802-827, 1995

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
%-------------------------------------------------------------------------------
function n=coloured_noise(TYPE,N,Fs)
if(nargin<1 || isempty(TYPE)) TYPE='WGN'; end 
if(nargin<2 || isempty(N)) N=256; end
if(nargin<3 || isempty(Fs)), Fs=50; end


switch TYPE
    % 1. White Gaussian Noise
  case 'WGN'
    w=randn(N,1);
    
    % 2. Coloured Gaussian Noise
    %    (1/f^alphas) noise process
  case 'CGN'
    alpha=5;
    w=noise_oneoverf(alpha,N,Fs,'FIR',1);
    
    % high-pass filter (0.3 Hz) to remove trend:
    w=filt_zerophase(w,Fs,1,0,0.3,201);
    w=w(1:N);
    
  otherwise
    error('Unknown noise-type.');
end

n=w;



function n=noise_oneoverf(alpha,N,fs,TYPE,Q_d,num_taps)
%-------------------------------------------------------------------------------
% Generate noise with a 1/(f^alpha) power spectra; code from [1]
%
% [1] Kasdin, N. Jeremy. "Discrete simulation of colored noise and stochastic processes
% and 1/f^α power law noise generation." Proc. IEEE, 83(5), 802-827, 1995
%-------------------------------------------------------------------------------
if(nargin<1 || isempty(alpha)) alpha=1; end
if(nargin<2 || isempty(N)) N=256; end
if(nargin<3 || isempty(fs)) fs=1; end
if(nargin<4 || isempty(TYPE)) TYPE='FREQ_DETERMIN'; end
if(nargin<5 || isempty(Q_d)) Q_d=1; end
if(nargin<6 || isempty(num_taps)) num_taps=N; end

DBplot=0;
DBwvd=0;
DBtfd=0;

TYPE=upper(TYPE);
white_gaus_noise=sqrt(Q_d).*randn(1,N);


% Code from Kasdin paper:
n_pts=N;
nn = N;
hfa=zeros(1,nn);    


%  /* find the deviation of the noise */
hfa(1)=1; 
for i=2:num_taps
    % /* generate the coefficients hk */ 
    hfa(i)=hfa(i-1)*(alpha/2+(i-2)) / (i-1) ;
end


%/* perform the discrete Fourier transform */ 
HFA=fft(hfa,nn);
WFA=fft(white_gaus_noise,nn);

n=ifft(HFA.*WFA);

n=real(n(1:N));

% Just for plotting:
if(DBplot) noise_spec=HFA.*WFA;  end


    





