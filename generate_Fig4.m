%-------------------------------------------------------------------------------
% generate_Fig4: Run the IF-estimation procedure to generate Fig. 4 from [1]
%
% Syntax: []=generate_Fig4()
%
% Inputs: 
%      - none
%
% Outputs: 
%     [] - none
%
% Example:
%     generate_Fig4;
%
% [1] B. Boashash, A. Ghazem, J.M. O' Toole, Time-frequency processing of non-stationary
% signals to aid diagnosis: highlights from medical applications, IEEE Signal Processing
% Magazine, 2013, in press

% John M. O' Toole, University College Cork
% Started: 10-05-2013
%-------------------------------------------------------------------------------
function []=generate_Fig4



%---------------------------------------------------------------------
% 0) parameters for plot
%---------------------------------------------------------------------
llimit=0; ulimit=2;
lower_freq_limit=0; upper_freq_limit=5;
FONT_TYPE='Times-Roman';
FONT_SIZE=14;
LINE_WIDTH=1.5;
Ntime=512;


% load the data:
b=load('PLED_example_epoch.mat');
x=b.x(1:1024); Fs=b.Fs;



%---------------------------------------------------------------------
% 1) generate the TFD of the signal
%---------------------------------------------------------------------
tf_orig=gen_TFD_EEG(x,Fs,Ntime);

figure(1); clf;
v_tfd90(tf_orig,x,Fs,lower_freq_limit,upper_freq_limit,FONT_TYPE,FONT_SIZE);
set(gca,'xlim',[0,4])
set(gca,'XTick',[0 1 2 3 4]);


%---------------------------------------------------------------------
% 2) Estimate and overlay the IF on top of the TFD
%---------------------------------------------------------------------
if_law=enhance_harmonics_estIF(x,Fs,Ntime);

hold on;
plot(if_law(:,2),if_law(:,1),'w','linewidth',3.5);
plot(if_law(:,2),if_law(:,1),'r','linewidth',2);




function v_tfd90(tfd,x,FS,f_bottom,f_top,FONT_TYPE,FONT_SIZE)
%---------------------------------------------------------------------
% PLOT TFD
%---------------------------------------------------------------------
[N,M]=size(tfd);
Ts=1/FS;
Ntime=length(x);

ntime=1:Ntime; ntime=ntime./FS;
n=linspace(ntime(1),ntime(end),N);
Mh=ceil(M/2);


tfd=tfd./max(tfd(:));
tfd=thres(tfd,0);


% this for spectrum:
kbottom=fix(f_bottom*Ntime*Ts);  ktop=fix(f_top*Ntime*Ts);
k=(kbottom+1:ktop);
N_k=length(kbottom:ktop);
% this for TFD:  (assming here that frequency domain is sampled with k/2NT)
kbottom=fix(f_bottom*2*M*Ts);  ktop=fix(f_top*2*M*Ts);
k_tfd=(kbottom+1:ktop);
f_tfd=k_tfd./(2*M*Ts);

imagesc(f_tfd,n,tfd(:,k_tfd)); axis xy;


ylim([0 20]);
% $$$ set(gca,'XTick',[0 4 8 12 16 20]);
set(gca,'YTick',[0 5 10 15 20]);
set(gca,'TickLength',[0.001, 0.025]);


ylabel('Time (secs)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
xlabel('Frequency (Hz)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
set(gca,'FontName',FONT_TYPE,'FontSize',FONT_SIZE);    




