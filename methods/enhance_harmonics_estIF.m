%-------------------------------------------------------------------------------
% enhance_harmonics_estIF: estimate IF using time-varying cepstrum approach
%
% Syntax: [if_law,if_law_samples,f_scale] = enhance_harmonics_estIF(x,Fs,Ntime,METHOD)
%
% Inputs: 
%     x      - time-domain signal (length N)
%     Fs     - sampling frequency
%     Ntime  - sampling frequency in time-direction of TFD
%     METHOD - either 'tfd' (proposed method) or 'tfdonly'
%
% Outputs: 
%     if_law         - estimate of IF
%     if_law_samples - estimate of IF in samples 
%     f_scale        - frequency scaling factor
%     tf             - TFD after homomorphic filtering
%     xp             - signal x after enhancing spikes (with NLEO)
%
% Example:
%     b=load('synth_signal_example_0dB.mat');
%     if_law=enhance_harmonics_estIF(b.x,b.Fs);
%
%     % plot:
%     figure(1); clf;  hold all;
%     plot(if_law(:,1),if_law(:,2)); 
%     plot( (1:length(b.true_IF))./b.Fs,b.true_IF);
%     legend('proposed','true IF');
%     xlim([10 30]); ylim([0 2]);
%     xlabel('time (seconds)'); 
%     ylabel('frequency (Hz)'); 
%

% John M. O' Toole, University College Cork
% Started: 09-05-2013
%-------------------------------------------------------------------------------
function [if_law,if_law_samples,f_scale,tf,xp]=enhance_harmonics_estIF(x,Fs,Ntime,METHOD)
if(nargin<3 || isempty(Ntime)) Ntime=[]; end
if(nargin<4 || isempty(METHOD)) METHOD='tfd'; end
if_law=[]; if_law_samples=[]; f_scale=[]; t_scale=[];



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




