%-------------------------------------------------------------------------------
% filter_zerophase: Simple implementation of zero-phase FIR filter using 'filtfilt'
%
% Syntax: x_filt=filter_zerophase(x,Fs,LP_fc,HP_fc,L_filt)
%
% Inputs: 
%     x,Fs,LP_fc,HP_fc,L_filt - 
%
% Outputs: 
%     x_filt - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 08-05-2013
%-------------------------------------------------------------------------------
function x_filt=filter_zerophase(x,Fs,LP_fc,HP_fc,L_filt)
if(nargin<2 || isempty(Fs)), Fs=1; end
if(nargin<3 || isempty(LP_fc)), LP_fc=[]; end
if(nargin<4 || isempty(HP_fc)), HP_fc=[]; end
if(nargin<5 || isempty(L_filt)), L_filt=[]; end

DBplot=0;

N=length(x);
if(isempty(L_filt))
    L_filt=make_odd(round(N/2-4)); 
end

%---------------------------------------------------------------------
% either bandpass, low-pass, or high-pass
%---------------------------------------------------------------------
if(~isempty(LP_fc) && ~isempty(HP_fc))
    filt_passband=[HP_fc/(Fs/2) LP_fc/(Fs/2)];
    b=fir1(L_filt,filt_passband);  
    
elseif(~isempty(LP_fc))
    filt_passband=LP_fc/(Fs/2);
    b=fir1(L_filt,filt_passband);

elseif(~isempty(HP_fc))
    L_filt=L_filt+1;
    filt_passband=HP_fc/(Fs/2);
    b=fir1(L_filt,filt_passband,'high');
    
else
    error('need to specify cut-off freuqency.');
end

  
% Do the filtering
x_filt=filtfilt(b,1,[x(:); zeros(L_filt,1)]);

x_filt=x_filt(1:N);

if(DBplot)
    figure(1); clf; 
    plot(x,'b'); hold on;
    plot(x_filt,'r');
end
