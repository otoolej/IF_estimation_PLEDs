%-------------------------------------------------------------------------------
% estIF_epoch_overlapadd: run method to extract IF from signal.  Assume that x is long and analysis
% needs to be done on an epoch by epoch basis. Overlap and add IF estimate.
%
% Syntax: [iflaw,t_scale,f_scale]=estIF_epoch_overlapadd(x,Fs,METHOD_TYPE)
%
% Inputs: 
%          x    - input signal (of length N)
%          Fs   - sampling frequency
%   METHOD_TYPE - 'tfd', 'tfdcosh', 'tfdonly', 'spike'
%
% Outputs: 
%     iflaw   - estimated IF law
%     t_scale - time scale factor
%     f_scale - frequency scale factor
%
% Example:
%     b=load('PLED_example_epoch.mat');
%     [iflaw,t_scale,f_scale]=estIF_epoch_overlapadd(b.x,b.Fs);
%
%     % plot:
%     figure(1); clf; 
%     n=1:length(iflaw);
%     plot(n.*t_scale,iflaw.*f_scale);
%     xlim([10 30]); ylim([0 5]);
%     xlabel('time (seconds)'); 
%     ylabel('frequency (Hz)'); 

% John M. O' Toole, University of Deusto
% Started: 21-02-2012
%-------------------------------------------------------------------------------
function [iflaw,t_scale,f_scale]=estIF_epoch_overlapadd(x,Fs,METHOD_TYPE)
if(nargin<2) error('Need input signal and sampling frequency.'); end
if(nargin<3 || isempty(METHOD_TYPE)) METHOD_TYPE='tfd'; end
% either 'tfd', 'tfdonly', or 'spike'

if(size(x,1)<2) x=x.'; end

iflaw=[]; 

DB=0;
DBverbose=1;
DBplot_iter=0;
DBplot_ONLY=0;


%---------------------------------------------------------------------
% Set parameters
%---------------------------------------------------------------------
EPOCH_LENGTH=20.48; % in seconds
% $$$ EPOCH_LENGTH=5.12; % (when using with spike-IF-estimation method) 
% $$$ EPOCH_LENGTH=40.96; % in seconds
OVERLAP=75; % in percent
WIN_TYPE='rect';
WIN_IF=1;
WIN_IF_TYPE='tukey'; WIN_IF_PARAM=0.2;

% if this fraction of the epoch is zero then don't process
FRAC_EPOCH_ZERO=0.85; 

%---------------------------------------------------------------------
% Assuming that N is too long to run TF analysis on.
%---------------------------------------------------------------------
if(mod(x,2))
    x=x(1:end-1);
end
if(length(x)<(EPOCH_LENGTH*Fs))
    disp('X shorter than epoch length; zero-padding...');
    x=zero_pad(x,ceil((EPOCH_LENGTH*Fs)));
    x=x.';
end

N=length(x);

% get epoch window and initialize variables:
[L_hop,L_epoch,win_epoch]=get_epoch_window(OVERLAP,EPOCH_LENGTH,WIN_TYPE,Fs);
N_epochs=floor( (N-L_epoch)/L_hop );
if(N_epochs<1) N_epochs=1; end



nw=0:L_epoch-1;
x_epoch=zeros(N_epochs,L_epoch);
y_test=zeros(N,1);
win_summed=zeros(N,1);

if(DB)
    dispVars(L_epoch,L_hop,N_epochs);
end

Ntime=L_epoch;
Ntime=512;

ratio_Ntime=L_epoch/Ntime;

if(mod(ratio_Ntime,1)~=0)
    error('need to have ratio_Ntime as integer.');
end
    

t_scale=ratio_Ntime/Fs;

nw_if=0:Ntime-1;

% need to scale for extracted IF, if TFD time-scale is different to signal:
Lhop_if=L_hop/ratio_Ntime;
N_if=floor(N/ratio_Ntime);


if_total=zeros(N_if,1);
if_weight=zeros(N_if,1);
if_win_weigth=zeros(N_if,1);



for k=1:N_epochs
    
    nf=mod(nw+(k-1)*L_hop,N);
    start_time=nf(1)/Fs;
    if(DBverbose)
        disp(['----- start time: ' num2str(start_time) ' ----']);
    end


    % epoch of the signal:
    x_epoch(k,:)=x(nf+1).*win_epoch;
    
    if(DBplot_iter)
        figure(23); clf;
        eeg_plot_simple(x_epoch(k,:),Fs,'epoch',23);
        drawnow;
    end

    
    %---------------------------------------------------------------------
    % extract the IF
    %---------------------------------------------------------------------
    if_law{k}=[]; if_law_samples=[];    
    ineg=find(x_epoch(k,:)==0);
    
    if(~DBplot_ONLY)

        % only process is signal is mostly non-zero:
        % (zeros are from the artefact removal process)
        if(length(ineg)<FRAC_EPOCH_ZERO*L_epoch)
            
            switch METHOD_TYPE
              case 'spike'
                [if_law{k},if_law_samples,f_scale]= ...
                    estIF_spike_method(x_epoch(k,:),Fs,Ntime);
                
              case {'tfd','tfdonly','tfdcosh'}
                [if_law{k},if_law_samples,f_scale]= ...
                    enhance_harmonics_estIF(x_epoch(k,:),Fs,Ntime,METHOD_TYPE);
            end
            
% $$$             if(DBplot_iter & ~isempty(tfd_filt))
% $$$                 plot_TFD_with_IFs(tfd_filt,x_epoch(k,:),Fs,if_law{k});
% $$$                 pause(1);
% $$$             end
        
        end
    

        nf_if=mod(nw_if+(k-1)*Lhop_if,N_if);
    
        %---------------------------------------------------------------------
        % overlap and add window on IF (adaptive to IF length)
        %---------------------------------------------------------------------
        if(WIN_IF)
            wif_pad=zeros(length(if_law_samples),1);
        
            in=find(if_law_samples>0);
            if(length(in)>0)
                wif=getWin(length(in),WIN_IF_TYPE,WIN_IF_PARAM);
                wif_pad(in)=wif;
            
                if(~isempty(if_law_samples))
                    if_law_samples=if_law_samples.*wif_pad;
                    if_win_weigth(nf_if+1)=if_win_weigth(nf_if+1)+wif_pad; % 
                end
            end

        end
    
        if(~isempty(if_law_samples))
            if_total(nf_if+1)=if_total(nf_if+1)+if_law_samples;
            
            ipos=find(if_law_samples>0);
            if_weight(nf_if(1)+ipos)=if_weight(nf_if(1)+ipos)+1;
        end
    
        
    
    

        y_test(nf+1)=y_test(nf+1)+x_epoch(k,:).';
        win_summed(nf+1)=win_summed(nf+1)+win_epoch;
        
    else
        %---------------------------------------------------------------------
        % Or just plot a TFD for the epoch:
        %---------------------------------------------------------------------
        
        tf=gen_TFD_EEG(x_epoch(k,:),Fs,Ntime,'lag-indep');
        figure(444); clf;
        wide_vtfd((tf),x_epoch(k,:),Fs,0,20);
    end
    
end
if(DBplot_ONLY)
    return;
end


y_test=y_test./win_summed;


if(WIN_IF)
    ipos=find(if_win_weigth>0);
    if_total(ipos)=if_total(ipos)./if_win_weigth(ipos);
else
    ipos=find(if_weight>0);
    if_total(ipos)=if_total(ipos)./if_weight(ipos);
end



%---------------------------------------------------------------------
% Post-processing. If artefacts were removed (and EEG set to zero), then 
% should negate the IF in this region too.
%---------------------------------------------------------------------
v=abs( sign(abs(x))-1 );

D=diff([0,v.',0]);
edge.start=find(D==1);
edge.last=find(D==-1)-1;

% set to zero, if zeros 'runs' > 1 second:
k=1; Li=length(edge.start);
for i=1:Li
    ni=edge.start(i):edge.last(i);
    if(length(ni)>Fs)
        if_total(ni)=0;
        k=k+1;
    end
end


iflaw=if_total;



%---------------------------------------------------------------------
% Debug: plot etc.
%---------------------------------------------------------------------
if(DB)
    figure(24); clf;
    eeg_plot_simple([x,y_test],Fs,{'orig','epoched'},24);
    
    figure(25); clf;
    n=(0:N_if-1).*t_scale;
    plot(n,if_total,n,if_weight,n,if_win_weigth);
end


if(DB)
    figure(27); clf;
    for k=1:N_epochs
        if(~isempty(if_law{k}))
            nf=mod(nw+(k-1)*L_hop,N);
            start_time=nf(1)/Fs;

            km=mod(k-1,7)+1;
            line_color=get(gca,'ColorOrder');
            plot(if_law{k}(:,1)+start_time,if_law{k}(:,2), ...
                 'Color',line_color(km,:));
            hold on;
        end
    end
end
    



    
