%-------------------------------------------------------------------------------
% gen_synth_signals_plus_noise: Generate many iterations of signal plus noise, 
% for given SNR
%
%
% Syntax: [sig_plus_noise,est_IFs]=gen_synth_signals_plus_noise(SNR, ...
%                                   N_iterations,N,Fs,NOISE_TYPE,WAVEFORM_MONOPHASIC);
%
% Inputs: 
%     SNR                 - signal-to-noise ratio (in dB)
%     N_iterations        - number of iterations
%     Fs                  - sample frequency
%     NOISE_TYPE          - noise type, either 'CGN' or 'WGN' (coloured/white 
%                           Gaussian noise)
%     WAVEFORM_MONOPHASIC - type of Duffing waveform, either 'monophasic' or
%                           'triphasic'
%
%
% Outputs: 
%     sig_plus_noise - signals
%     est_IFs        - instantaneous frequency of signal
%
% Example:
%     snr=0; Fs=50; 
%     N_iters=10; WAVEFORM_MONOPHASIC=0;
%     [sigs,IFs]=gen_synth_signals_plus_noise(snr,N_iters,Fs,'CGN',WAVEFORM_MONOPHASIC);
%
%     figure(1); clf; Fs=50; 
%     n=(1:size(sigs,2))./Fs; 
%     subplot(2,1,1); plot(n,sigs(1,:)); xlim([5 35]);
%     xlabel('time (seconds)');  ylabel('amplitude');
%     subplot(2,1,2); plot(n,IFs(1,:)); xlim([5 35]);
%     xlabel('time (seconds)');  ylabel('frequency (Hz)');
%

% John M. O' Toole, University College Cork
% Started: 08-04-2013
%-------------------------------------------------------------------------------
function [sig_plus_noise,est_IFs]=gen_synth_signals_plus_noise(SNR,N_iterations, ...
                                                  Fs,NOISE_TYPE,WAVEFORM_MONOPHASIC)
if(nargin<1 || isempty(SNR)), SNR=[]; end
if(nargin<2 || isempty(N_iterations)), N_iterations=20; end
if(nargin<3 || isempty(Fs)), Fs=50; end
if(nargin<4 || isempty(NOISE_TYPE)), NOISE_TYPE='CGN'; end
if(nargin<5 || isempty(WAVEFORM_MONOPHASIC)), WAVEFORM_MONOPHASIC=1; end


DBplot=0;
DBplot_TFD=0;

t_parameters;

%---------------------------------------------------------------------
% 1) generate the different 'discharge' waveforms 
%---------------------------------------------------------------------
t=load(DUFFING_WAVEFORMS_FILE);
N=2049;

if(WAVEFORM_MONOPHASIC)
    duff=t.duff_monophasic;
else
    duff=t.duff_triphasic;
end


N_duff=size(duff,1);


sig_plus_noise=zeros(N_iterations*N_duff,N);
sig_no_noise=zeros(N_iterations*N_duff,N);
est_IFs=zeros(N_iterations*N_duff,N);

% can load signals (noise free) from file:
sig_st_mat=[];
    

nn=1;
for p=1:N_duff
    
    dispVars(p);
    for n=1:N_iterations
        
        %---------------------------------------------------------------------
        % get train of duffing waveforms
        %---------------------------------------------------------------------
        [sp,imp_train,est_IF]=duffing_train(duff(p,:),N,Fs);
        est_IFs(nn,:)=est_IF;
        sp=scale_s(sp,1);


        %---------------------------------------------------------------------
        % get noise
        %---------------------------------------------------------------------
        if(~isempty(SNR))
            no=coloured_noise(NOISE_TYPE,N,Fs);
            no=scale_s(no,1);
        
            no=scale_n(no,getEnorm(sp),SNR);
            
            sig_plus_noise(nn,:)=sp+no(:);
            sig_plus_noise(nn,:)=scale_s(sig_plus_noise(nn,:),1);
        else
            sig_no_noise(nn,:)=sp;
        end
        
        
        nn=nn+1;
    end
end

if(~isempty(sig_st_mat))
    est_IFs=sig_st_mat.est_IFs; 
end


%---------------------------------------------------------------------
% PLOT
%---------------------------------------------------------------------
if(DBplot)
  figure(18); clf;
  plot( (1:N)./Fs, sp);
end

if(DBplot_TFD)
   dt=dtfd_sep1(sp,{11,'hamm',0,1},{1024,'bartlett'},floor(N/4),N*2); 
   figure(19); clf; 
   wide_vtfd(dt,sp,Fs,0,10);
end



function s=scale_s(s,Enew)
%---------------------------------------------------------------------
% Scale signal so that that energy of signal is equal to Enew 
%---------------------------------------------------------------------
Es=getEnorm(s);
s=s.*sqrt(Enew/Es);


function E=getEnorm(s)
%---------------------------------------------------------------------
% energy of signal
%---------------------------------------------------------------------
E=sum( abs(s).^2 )/length(s);


function new_noise=scale_n(noise,E_signal,SNR)
%---------------------------------------------------------------------
% Scale noise to required SNR level (in dB)
%---------------------------------------------------------------------
N=length(noise);
E_noise=getEnorm(noise);

E_new_noise = E_signal/(10^(SNR/10));

new_noise=noise.*sqrt(E_new_noise/E_noise);  

