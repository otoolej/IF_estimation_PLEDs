%--------------------------------------------------------------------------------
% General function to calculate a window function. (update of calWin.m)
%
% function win = getWin( win_length, win_type, aux_param, freq_domain_rep )
%
%
% Returns vector with negative then positive indices of window
%
%
% STARTED: 05-01-2008
%--------------------------------------------------------------------------------

function win = getWin( win_length, win_type, aux_param, freq_domain_rep )
if( nargin<3 ) aux_param=[]; end
if( nargin<4 ) freq_domain_rep=0; end
DB=0;



switch win_type
 case 'delta'
  win=zeros(win_length,1);
  wh = floor(win_length/2);
  win(wh+1)=1;
 case 'rect'
  win(1:win_length) = 1;
 case 'bart'
  win = bartlett( win_length );
 case {'hamm', 'hamming'}
  win = hamming( win_length );
 case { 'hann', 'hanning' }
  win = hanning( win_length );
 case 'tukey'
  %% Problem with Octave's tukeywin.m for N odd,
  %% use Matlab's one instead.
  if(nargin<3 || isempty(aux_param))      
    win = tukeywin( win_length );
  else
    win = tukeywin( win_length, aux_param );
  end

 case 'gauss'
  if(nargin<3 || isempty(aux_param))      
    win = gausswin( win_length );
  else
    win = gausswin( win_length, aux_param );
  end
 case 'cosh'
  win_hlf = fix( win_length / 2);

  if(nargin<3 || isempty(aux_param))
    aux_param = 0.01;
  end
  for m = -win_hlf:win_hlf
    win(mod(m,win_length)+1) = cosh( m ).^( -2 * aux_param );
  end
  win = fftshift(win);

  %% Nonsymmetry TEST window functions.
 case 'randn'
  win=randn(win_length,1);
 case 'sin'
  f_0=4.44;
  win=sin( 2*pi*f_0.*(0:win_length-1)./win_length );
 case 'rload'
  % LOAD this one from a data file:
  % (vector of 1024 randn sample points stored)
  if(win_length>1024)
    error('Max. length for this window 1024');
  end
  % $$$   b=load('WRANDN.dat','wrand');
  FILE='~/work/software/octave_matlab/octave_only/WRANDN.dat';
  b=load(FILE,'wrand');
  win=b.wrand(1:win_length);


  %% Conjugate symmetric random test window
 case 'randnsysm'
  wh=floor(win_length/2);
  n=0:wh;
  win(n+1)=randn(wh+1,1)+j.*randn(wh+1,1);

  %% need to force this if want DFT to be real valued
  win=force_conj_symm(win,win_length);


  %% And this is the way we have been returning the windows:
  win_shift=circshift(win(:),floor(win_length/2));
  win=win_shift;


 case 'sinsysm'
  wh=floor(win_length/2);
  n=0:wh;

  f_0=4.44; f_1=3.7;
  win=sin( 2*pi*f_0.*n./win_length ) + j.*sin( 2*pi*f_1.*n./win_length );



  %% need to force this if want DFT to be real valued
  win=force_conj_symm(win,win_length);

  %% And this is the way we have been returning the windows:
  win_shift=circshift(win(:),floor(win_length/2));
  win=win_shift;



 otherwise
  error(['Unknown window type ' win_type]);
end



%---------------------------------------------------------------------
% When want the DFT of win to be equal to the 'usual' window function.
%---------------------------------------------------------------------
if( freq_domain_rep==1 )

  if(DB) dispVars('Returning DFT of window'); end
  win_hlf=ceil(win_length/2);


  win=circshift(win(:),ceil(win_length/2));
  win=fft(win);

  % Then shift in frequency domain
  % (which is the antithesis of above shift)
  win=circshift(win(:),floor(win_length/2));


  %% modulate to componsate for the circular-shift.
  %% OK. Should fix up at some stage.
% $$$   l=0:win_length-1;
% $$$   mo=exp(-j*2*pi*floor(win_length/2).*l./win_length).';
% $$$   win=win(:).*mo;
  %% This is the long way:
% $$$   win=ifft( circshift(fft(win)(:),floor(win_length/2)) );



end






function win=force_conj_symm(win,win_length)
DB=1;
wh=floor(win_length/2);
even_factor=1-rem(win_length,2);

n=1:wh-even_factor;
win(win_length-n+1)=conj(win(n+1));

win(1)=real(win(1));
if(even_factor==1)
  win(wh+1)=real(win(wh+1));
end


if(DB)
  if( isreal_fn(fft(win))==0 )
    error('Not conjugate symmetrical.');
  end
end
