%---------------------------------------------------------------------
% estimate IFs (3 methods) at different SNR values:
%---------------------------------------------------------------------

% range of signal-to-noise ratios:
snrs=-20:2.5:20;
% 20 iterations for each of the 5 Duffing waveforms:
N_iterations=20; 


% generate mean-squared error for each of the two synthetic datasets:
WAVEFORM_MONOPHASIC=0;
for n=1:length(snrs)
    err_mono(n)=if_estimation_err(snrs(n),N_iterations,WAVEFORM_MONOPHASIC);
end
plot_results_SNR_ifests(err_mono,WAVEFORM_MONOPHASIC);

WAVEFORM_MONOPHASIC=1;
for n=1:length(snrs)
    err_tri(n)=if_estimation_err(snrs(n),N_iterations,WAVEFORM_MONOPHASIC);
end
plot_results_SNR_ifests(err_tri,WAVEFORM_MONOPHASIC);
