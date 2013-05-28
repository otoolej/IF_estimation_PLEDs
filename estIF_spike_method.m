%-------------------------------------------------------------------------------
% estIF_spike_method: Estimate the IF using a spike-enhancing method; idea from [1]
%
% Syntax: [if_law,if_law_samples,f_scale]=estIF_spike_method(x,Fs,Ntime)
%
% Inputs: 
%     x         - time-domain signal (length N)
%     Fs        - sampling frequency
%     Ntime     - downsampling the IF from N to Ntime 
%
% Outputs: 
%     if_law         - estimate of IF 
%     if_law_samples - estimate of IF in samples 
%     f_scale        - frequency scaling factor
%
% Example:
%     b=load('PLED_example_epoch.mat');
%     [iflaw,~,f_scale]=estIF_spike_method(b.x,b.Fs);
%
%     % plot:
%     figure(1); clf; 
%     plot(iflaw(:,1),iflaw(:,2));
%     xlim([10 30]); ylim([0 5]);
%     xlabel('time (seconds)'); 
%     ylabel('frequency (Hz)'); 
%
%
% [1] Deburchgraeve, W., Cherian, P. J., De Vos, M., Swarte, R. M., Blok, J. H., Visser,
% G. H., Govaert, P., et al. (2008). Automated neonatal seizure detection mimicking a
% human observer reading EEG. Clinical Neurophysiology, 119(11), 2447–2454.


% John M. O' Toole, University of Deusto
% Started: 25-01-2012
%-------------------------------------------------------------------------------
function [if_law,if_law_samples,f_scale]=estIF_spike_method(x,Fs,Ntime)
if(nargin<3 || isempty(Ntime)) Ntime=length(x); end

DBplot=0;


if_law=[]; if_law_samples=[]; tf=[];

f_scale=Fs;
N=length(x);
t_scale=((N/Fs)/Ntime);


RE_TFD=0;

ULIMIT_FREQ_EST=5; % in Hz (remove any estimate above this)


% a) enhance the spikes (assuming spike-like signals)
% $$$ [xp,x_peaks]=spike_emphasis(x,Fs,'DEBURCH');
[xp,x_peaks]=locate_spikes(x,Fs,'DEBURCH',1);

nzeros=find(x_peaks>min(x_peaks));
xnz=zeros(N,1);
xnz(nzeros)=x(nzeros);

% b) interpolate between spikes to get IF
ipeaks=find(abs(xnz)>0);
rr=diff(abs(ipeaks));

if(DBplot) rro=rr; ipeakso=ipeaks; end


% trim short periods (below 1/5Hz)
ishort=find(rr<(1/ULIMIT_FREQ_EST)*Fs);
rr(ishort)=[];
ipeaks(ishort)=[];

ff=1./rr;

if(DBplot) ffo=1./rro; end


ipeaks_p=[];
if(length(ff)>1)
    in=ipeaks(2):ipeaks(end);
    ipeaks_p=ipeaks(2:end);


    ff_interp=interp1( ipeaks_p, ff, in, 'cubic' );    


    ff_interp_all=zeros(1,N);
    ff_interp_all(in)=ff_interp;


    if_law_samples=ff_interp_all;

    % downsample to Ntime:
    if_law_samples=if_law_samples(1:(N/Ntime):end).';

    if_law(:,1)=find(if_law_samples>0);
    if_law(:,2)=if_law_samples(if_law(:,1));

    if_law(:,1)=if_law(:,1).*t_scale;
    if_law(:,2)=if_law(:,2).*f_scale;
end








