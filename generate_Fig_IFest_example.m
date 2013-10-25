%-------------------------------------------------------------------------------
% generate_Fig_IFest_example: Plot TFD and IF estimates for PLED example
%
% Syntax: []=generate_Fig_IFest_example(fig)
%
% Inputs: 
%     fig, - either: 'A','B','C','D', or 'E'
%
% Outputs: 
%     [] - none
%
% Example:
%     generate_Fig_IFest_example('B');
%

% John M. O' Toole, University College Cork
% Started: 26-06-2013
%-------------------------------------------------------------------------------
function []=generate_Fig_IFest_example(fig)
if(nargin<1 || isempty(fig)) fig='a'; end
if(nargin<2 || isempty(x)) x=[]; end
if(nargin<3 || isempty(Fs)) Fs=100; end

fig=lower(fig);

% load parameters:
t_parameters;

% examine IF in this range:
llimit=0; ulimit=1;
lower_freq_limit=0; upper_freq_limit=3;


% load PLED signal from .mat file
b=load([DATA_DIR 'PLED_example_20sec_epoch.mat']);
x=b.x; Fs=50;


%---------------------------------------------------------------------
% time-domain signal
%---------------------------------------------------------------------
if( strcmp(fig,'a')==1 )
  p_EEG(9,x,Fs,FONT_TYPE,FONT_SIZE);
end


%---------------------------------------------------------------------
% TFD with raw EEG signal
%---------------------------------------------------------------------
if( strcmp(fig,'b')==1 )
  tf_orig=g_TFD(x,Fs);

  % for plotting normalise and remove anything below 0:
  tf=tf_orig./max(max(tf_orig));
  tf(tf<0)=0;

  figure(11); clf;
  v_tfd(tf,x,Fs,lower_freq_limit,upper_freq_limit,FONT_TYPE,FONT_SIZE);
end


%---------------------------------------------------------------------
% TFD with spike signal (separable kernel)
%---------------------------------------------------------------------
if( strcmp(fig,'c')==1 )
  tf=loglag_filter_TFD(x,Fs);    

  figure(16); clf;    
  v_tfd(tf,x,Fs,lower_freq_limit,upper_freq_limit,FONT_TYPE,FONT_SIZE);
end


%---------------------------------------------------------------------
% Estimate IF from TFD
%---------------------------------------------------------------------
if( strcmp(fig,'d')==1 )

  [tf,t_scale,f_scale]=g_TFD(x,Fs);  
  tf=tf./max(max(tf));

  if_law=extract_IFtracks(tf,length(x),Fs);

  p_IF_tracks(21,if_law,t_scale,f_scale,...
              llimit,ulimit,FONT_TYPE,FONT_SIZE);
end

%---------------------------------------------------------------------
% Estimate IF from TFD and cep-TFD
%---------------------------------------------------------------------
if( strcmp(fig,'e')==1 )
  tf=loglag_filter_TFD(x,Fs);    

  [if_law,~,f_scale,t_scale]=extract_IFtracks(tf,length(x),Fs);


  p_IF_tracks(22,if_law,t_scale,f_scale,...
              llimit,ulimit,FONT_TYPE,FONT_SIZE);
end



function tf=loglag_filter_TFD(x,Fs)
%---------------------------------------------------------------------
%% generate log-lag filtered TFD
%---------------------------------------------------------------------
xp=locate_spikes(x,Fs,'NEO',1);
xp=xp-mean(xp);
tf_xp=g_TFD(xp,Fs);

tf_filt=loglag_TFfilter(tf_xp,Fs,0.00025);        

tf=tf_filt; x=xp;        
tf=tf-mean(tf(:));
tf=tf./max(max(tf));
tf(tf<0)=0;



function [dt,t_scale,f_scale]=g_TFD(x,Fs)
%---------------------------------------------------------------------
%% TFD from here as parameters may change in gen_TFD_EEG.m
%---------------------------------------------------------------------
Ntime=512; 
Nfreq=8192*4;    

dur=length(x)/Fs;
t_scale=(dur/Ntime);
f_scale=(1/Nfreq)*(Fs/2);

N=length(x);

win_lag={get_odd_int(15*Fs)-1,'hamm'}; 


f_dopp_Hz=1.15; % length of doppler window in Hertz    

f_dopp_Hz=0.7; % length of doppler window in Hertz    
L_dopp=make_odd( floor(f_dopp_Hz*(N/Fs)) );
win_dopp={L_dopp,'hann',0.01,1};


dt=dtfd_sep2(x,win_dopp,win_lag,Ntime,Nfreq);



function v_tfd(tfd,x,FS,f_bottom,f_top,FONT_TYPE,FONT_SIZE)
%---------------------------------------------------------------------
%% PLOT TFD
%---------------------------------------------------------------------
[N,M]=size(tfd);
Ts=1/FS;
Ntime=length(x);


ntime=1:Ntime; ntime=ntime./FS;
n=linspace(ntime(1),ntime(end),N);
Mh=ceil(M/2);


% this for spectrum:
kbottom=fix(f_bottom*Ntime*Ts);  ktop=fix(f_top*Ntime*Ts);
k=(kbottom+1:ktop);
N_k=length(kbottom:ktop);
% this for TFD:  (assming here that frequency domain is sampled with k/2NT)
kbottom=fix(f_bottom*2*M*Ts);  ktop=fix(f_top*2*M*Ts);
k_tfd=(kbottom+1:ktop);
f_tfd=k_tfd./(2*M*Ts);


imagesc(n,f_tfd,tfd(:,k_tfd).'); axis xy;


xlim([0 20]);
set(gca,'XTick',[0 20]);
set(gca,'YTick',[0:3]);

set(gca,'TickLength',[0.001, 0.025]);

% add lines where we expect the IF:
LINE_WIDTH=2;

line([0 20],[0.2 0.2],'color','y','LineWidth',LINE_WIDTH, ...
     'linestyle','-.');
line([0 20],[0.5 0.5],'color','y','LineWidth',LINE_WIDTH, ...
     'linestyle','-.');


xlabel_font('TIME (seconds)',FONT_TYPE,FONT_SIZE);
ylabel_font('FREQUENCY (Hz)',FONT_TYPE,FONT_SIZE);
set(gca,'FontName',FONT_TYPE,'FontSize',FONT_SIZE);    




function p_EEG(fig_num,x,Fs,FONT_TYPE,FONT_SIZE)
%---------------------------------------------------------------------
%% PLOT IFs
%---------------------------------------------------------------------
h=figure(fig_num); 

FONT_SIZE=FONT_SIZE+8;

scrsz = get(0,'ScreenSize');
set(h,'Position',[1 scrsz(4)/1 scrsz(3)*0.8 scrsz(4)/4]); 
clf; 

LINE_WIDTH=1.4;


n=0:length(x)-1; t=n./Fs;
hp=plot(t,x);

set(hp,'Color', [0 0 .5],'LineWidth',LINE_WIDTH);

xlim([0 20]);
ylim([-300 100]);

% hack: draw the grid lines on manually:
y_limits=ylim; 
cc=[1 1 1]*0.7;
grid_width=1.2;

for i=2:2:21
    L=line([i i],[y_limits(1) y_limits(2)],'color',cc, ...
         'linewidth',grid_width);
    uistack(L, 'bottom');
end

y_lims=ylim;
line([2 4], [y_lims(1)+4 y_lims(1)+4],'color','k','linewidth',4);
line([2 2], [y_lims(1)+2 y_lims(1)+102],'color','k','linewidth',4);
text(2.3,-290,'2 seconds','FontName',FONT_TYPE,'FontSize',FONT_SIZE,...
     'VerticalAlignment','bottom','HorizontalAlignment','left');
text(1.9,-290,'100 \muV','FontName',FONT_TYPE,'FontSize',FONT_SIZE,...
     'VerticalAlignment','bottom','HorizontalAlignment','left',...
     'rotation',90);

set(gca,'box','off');
set(gca,'xtick',[]); set(gca,'ytick',[]);
axis off;

set(gca,'FontName',FONT_TYPE,'FontSize',FONT_SIZE);    







function p_IF_tracks(fig_num,if_law,t_scale,f_scale,llimit,...
                     ulimit,FONT_TYPE,FONT_SIZE)
%---------------------------------------------------------------------
%% PLOT IFs
%---------------------------------------------------------------------
h=figure(fig_num); 

scrsz = get(0,'ScreenSize');
set(h,'Position',[1 scrsz(4)/1 scrsz(3)*0.38 scrsz(4)/4]); 
clf; 

LINE_WIDTH=1.5;

hp=plot(if_law(:,1),if_law(:,2),'.');

ylim([llimit,ulimit]);
set(gca,'XTick',[0 20]);
set(gca,'YTick',[0 0.5 1]);

xlim([0 20]);

xlabel_font('TIME (seconds)',FONT_TYPE,FONT_SIZE);
ylabel_font('FREQUENCY (Hz)',FONT_TYPE,FONT_SIZE);
set(gca,'FontName',FONT_TYPE,'FontSize',FONT_SIZE);    

set(hp,'Color', [0 0 .5],'LineWidth',LINE_WIDTH);
set(gca,'TickDir','out');


set(gca,'box','off');
% hack: draw the grid lines on manually:
x_limits=xlim; 
cc=[1 1 1]*0.7;
grid_width=2;

L=rectangle('position',[0 0.2 20 0.3],'EdgeColor','none','FaceColor',cc);
uistack(L, 'bottom');




