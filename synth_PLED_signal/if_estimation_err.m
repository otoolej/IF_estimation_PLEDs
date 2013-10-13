%-------------------------------------------------------------------------------
% if_estimation_err: Estimate and compare the IF from different three methods. Test
% signals are read in from .mat files
%
% WARNING: need at least 2GB RAM (PC memory)
%
%
% Syntax: [err_st]=if_estimation_err(SNR,N_iterations,WAVEFORM_MONOPHASIC)
%
% Inputs: 
%     SNR                 - signal-to-noise ratio (in dB)
%     N_iterations        - number of iterations
%     WAVEFORM_MONOPHASIC - type of Duffing waveform, either 'monophasic' or
%                           'triphasic'
%
% Outputs: 
%     err_st - structure with mean-squared error for the 3 methods
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 21-05-2013
%-------------------------------------------------------------------------------
function [err_st]=if_estimation_err(SNR,N_iterations,WAVEFORM_MONOPHASIC)
if(nargin<1 || isempty(SNR)), SNR=[]; end
if(nargin<2 || isempty(N_iterations)), N_iterations=100; end
if(nargin<3 || isempty(WAVEFORM_MONOPHASIC)), WAVEFORM_MONOPHASIC=0; end


DBplot_cont=0;
DBplot=1;
DBverbose=1;
SAVE_RESULTS=0;


t_parameters;

%---------------------------------------------------------------------
% 1. generate the synthetic data:
%---------------------------------------------------------------------
Fs=50; 
[sig_plus_noise,est_IFs]=gen_synth_signals_plus_noise( ...
    SNR,N_iterations,Fs,'CGN',WAVEFORM_MONOPHASIC);
    

L=size(sig_plus_noise,1);
err_tfd=zeros(L,1); err_tfdonly=zeros(L,1); err_spike=zeros(L,1);
if(DBverbose)
    fprintf(' -- SNR = %d -- \n',SNR); 
    fprintf('  | proposed | TFD only | spike |\n');     
    fprintf('---------------------------------\n');         
end


for ii=1:L
    if(DBverbose), fprintf('%d',ii); end

    x=sig_plus_noise(ii,1:end-1); 
    
    IF_EST_METHOD='tfd';
    [iflaw,t_scale,f_scale]=estIF_epoch_overlapadd(x,Fs,IF_EST_METHOD);
    iflaw_tfd=iflaw.*f_scale;

    IF_EST_METHOD='tfdonly';
    [iflaw,t_scale,f_scale]=estIF_epoch_overlapadd(x,Fs,IF_EST_METHOD);
    iflaw_tfdonly=iflaw.*f_scale;

    IF_EST_METHOD='spike';
    [iflaw,t_scale,f_scale]=estIF_epoch_overlapadd(x,Fs,IF_EST_METHOD);
    iflaw_spike=iflaw.*f_scale;


    %---------------------------------------------------------------------
    % interpolate IF laws so the same as 'true' IF (if using time decimation in DTFD)
    %---------------------------------------------------------------------
    n=(1:length(iflaw_tfd)).*t_scale; in=(1:length(est_IFs))./Fs;
    iflaw_tfd=interp1(n,iflaw_tfd,in,'pchip');
    iflaw_tfdonly=interp1(n,iflaw_tfdonly,in,'pchip');
    iflaw_spike=interp1(n,iflaw_spike,in,'pchip');


    %---------------------------------------------------------------------
    % compare IF: mean squared error
    %---------------------------------------------------------------------

    % measure only the middle 20 seconds:
    n=ceil(10*Fs):floor(30*Fs);
    ns=ceil(10*Fs):floor(30*Fs);


    err_tfd(ii)=mean( abs(est_IFs(ii,n)-iflaw_tfd(n)).^2 ); %./norm_y1;
    err_tfdonly(ii)=mean( abs(est_IFs(ii,n)-iflaw_tfdonly(n)).^2 ); %./norm_y1;
    err_spike(ii)=mean( abs(est_IFs(ii,n)-iflaw_spike(n)).^2 ); %./norm_y1;

    if(DBverbose)
        fprintf(' | %1.5f | %1.5f | %1.5f |\n', ... 
                err_tfd(ii),err_tfdonly(ii),err_spike(ii));
    end
    
    
    if(DBplot_cont)
        figure(6); clf; hold all;
        t=(1:length(iflaw_tfd))./Fs;
        plot(t,iflaw_tfd,'g'); 
        plot(t,iflaw_tfdonly,'b');   
        plot(t,iflaw_spike,'r');     
        plot(t,est_IFs(ii,:), 'k' );
        legend('proposed','tfd','spike','true IF');
        
        drawnow;
    end
    
end

if(DBverbose), fprintf('\n'); end

err_st.SNR=SNR;
err_st.method_tfd=err_tfd;
err_st.method_tfdonly=err_tfdonly;
err_st.method_spike=err_spike;


if(SAVE_RESULTS)
    if(~WAVEFORM_MONOPHASIC), IMP_STR=[]; end
    
    fname_out=[RES_DIR IMP_STR 'if_err_' num2str(SNR) '_dB_' ...
               datestr(now,'dd-mmm-yyyy_HH:MM:SS') '.mat'];
    save(fname_out,'err_st');
end

  
