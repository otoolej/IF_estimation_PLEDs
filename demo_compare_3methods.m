%-------------------------------------------------------------------------------
% demo_compare_3methods: Compare 3 methods for estimate the instantaneous frequency
% (IF) of a test signal
%
% Syntax: []=demo_compare_3methods()
%
% Inputs: 
%      - 
%
% Outputs: 
%      - 
%
% Example:
%     demo_compare_3methods;
%

% John M. O' Toole, University College Cork
% Started: 26-10-2013
%-------------------------------------------------------------------------------
function []=demo_compare_3methods()


% 0. load the test (synthetic) signals:
b=load('synth_signal_example_0dB.mat');
x=b.x; Fs=b.Fs;

% 1. try the 3 different methods:
[iflaw_spike,t_scale_spike,f_scale_spike]=estIF_epoch_overlapadd(x,Fs,'spike');
[iflaw_tfdonly,t_scale,f_scale]=estIF_epoch_overlapadd(x,Fs,'tfdonly');
[iflaw_tfd,t_scale,f_scale]=estIF_epoch_overlapadd(x,Fs,'tfd');

% 2. plot:
figure(2); clf;  hold all;
t=(1:length(iflaw)).*t_scale;
plot(t,iflaw_tfd.*f_scale);
plot(t,iflaw_tfdonly.*f_scale);
plot((1:length(iflaw)).*t_scale_spike,iflaw_spike.*f_scale_spike);
n=1:length(true_IF);
plot(n./Fs,true_IF);
legend('proposed','TFD only','spike detection','true IF');
xlim([10 30]); ylim([0 3]);
xlabel('time (seconds)'); 
ylabel('frequency (Hz)'); 
