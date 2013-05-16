%-----------------------------------------------------------------------------------------
% FFT ROUTINE for input data with conjugate symmetry.
% 
% output should be real-valued.
%-----------------------------------------------------------------------------------------
function Y=fft_conj_symm(x)
Y=fft(x);
