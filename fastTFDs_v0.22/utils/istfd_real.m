%-------------------------------------------------------------------------
% Check if DTFD is real valued.
%
%-------------------------------------------------------------------------
function tf=istfd_real(tf,s1,fromfile)

E=sum(abs(s1).^2) * 1e-12;
if( max( abs(imag(tf(:))) ) > E )
  if(nargin>2)
    disp(['**-&- WARNING: TFD from ' fromfile ' not real -&-**']);
  else
    disp('**-&- WARNING: TFD not real -&-**');
  end
else
  tf=real(tf);
end
