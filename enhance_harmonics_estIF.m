%-------------------------------------------------------------------------------
% enhance_harmonics_estIF: estimate IF using time-varying cepstrum approach
%
% Syntax: if_law = enhance_harmonics_estIF(x,Fs,Ntime,METHOD)
%
% Inputs: 
%     x      - time-domain signal (length N)
%     Fs     - sampling frequency
%     Ntime  - sampling frequency in time-direction of TFD
%     METHOD - either 'tfd' (proposed method) or 'tfdonly'
%
% Outputs: 
%     if_law - IF law
%
% Example:
%    b=load('PLED_example_epoch.mat');
%    if_law=enhance_harmonics_estIF(b.x,b.Fs);
%
%    figure(1); clf; 
%    plot(if_law(:,1),if_law(:,2)); ylim([0 5]);
%    xlabel('Time (secs)'); ylabel('Frequency (Hz)');
%

% John M. O' Toole, University College Cork
% Started: 09-05-2013
%-------------------------------------------------------------------------------
function if_law=enhance_harmonics_estIF(x,Fs,Ntime,METHOD)
if(nargin<3 || isempty(Ntime)) Ntime=[]; end
if(nargin<4 || isempty(METHOD)) METHOD='tfd'; end


if_law=[]; if_law_samples=[]; f_scale=[]; t_scale=[];
DBplot=0;


IF_EST_METHOD='MQ';
LOGLAG_FILTER=0.00025;

if(strcmp(METHOD,'tfdonly')==1)
    LAG_METHOD='none';
else
    LAG_METHOD='loglag';
end



%---------------------------------------------------------------------
% a) enhance the spikes (assuming spike-like signals)
%---------------------------------------------------------------------
if(strcmp(METHOD,'tfdonly')==0)
    xp=locate_spikes(x,Fs,'NEO',1);
else
    xp=x;
end


%---------------------------------------------------------------------
% b) generate TFD
%---------------------------------------------------------------------
xp=xp-mean(xp);
if(strcmp(LAG_METHOD,'none'))
    xp=x;
end
tf=gen_TFD_EEG(xp,Fs,Ntime);


if(DBplot)
    figure(43); clf;
    wide_vtfd(thres(tf,0),x,Fs,0,10);
end
    


%---------------------------------------------------------------------
% c) filter in the time-varying 'cepstrum' domain
%---------------------------------------------------------------------
if( strcmp(LAG_METHOD,'loglag')==1 )
    tf=loglag_TFfilter(tf,Fs,LOGLAG_FILTER);        
end


%---------------------------------------------------------------------
% d) extract IF from tfd:
%---------------------------------------------------------------------
[if_law,if_law_samples,f_scale,t_scale]= ...
    extract_IFtracks(tf,length(x),Fs);


if(DBplot)
    figure(44); clf;
    wide_vtfd(thres(tf,0),x,Fs,0.2,[]);
end



