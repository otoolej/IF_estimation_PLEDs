%-------------------------------------------------------------------------------
% Get discrete analytic signal.  See [1] for details.
%
% USE: z=get_analytic(s)
%
% INPUT:
%      s = real-valued signal
%
% OUTPUT:
%      z = complex-valued analytic signal
%
% EXAMPLE
%      N=64;  f=0.1; N2=2*N;
%      s=cos( 2*pi*f.*(0:N-1));     
%      z=get_analytic(s);
%      plot(1:N2,real(z),1:N2,imag(z));
%      legend('real','imaginary');
%
% 
% [1] J. M. O' Toole, M. Mesbah, and B. Boashash, "A New Discrete Analytic Signal 
%     for Reducing Aliasing in the Discrete Wigner-Ville Distribution", IEEE Trans. 
%     on Signal Processing,  vol. 56, no. 11, pp. 5427-5434, Nov. 2008.


% Copyright (C) 2007,2008 John M. O' Toole, The University of Queensland
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%-------------------------------------------------------------------------------
function z=get_analytic(s1);
s1=real(s1);
N=length(s1); N2=2*N;

% 1. zero-pad N-point signal to 2N
s1=[s1(:); zeros(N,1)];
S1=fft(s1);

% 2. Get analytic signal of 2N-point signal
H=zeros(N2,1); 
H([1 N+1])=1;
H(2:N)=2;
z_cb=ifft(S1.*H);


% 3. Force the second half of the time-domain signal to zero.
z=[z_cb(1:N); zeros(N,1)];







