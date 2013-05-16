%-----------------------------------------------------------------------------------------
% IFFT for input data with conjugate symmetry.
% 
% output is should be real-valued.
%-----------------------------------------------------------------------------------------
function Y=ifft_conj_symm(x)
Y=ifft(x);
