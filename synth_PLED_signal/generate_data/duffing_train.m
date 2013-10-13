%-------------------------------------------------------------------------------
% duffing_train: Generate nonstionary sequeuence from one Duffing waveform
%
% Syntax: [sp,imp_train,est_IF]=duffing_train(duff,N,Fs)
%
% Inputs: 
%     duff - waveform (output of Duffing oscillator)
%     N    - signal length
%     Fs   - sample frequency
%
% Outputs: 
%     sp         - signal, a nonstationary sequence of Duffing waveforms
%     imp_train  - nonstationary impulse train
%     est_IF     - instantaneous frequency of signal
%
% Example:
%     t_parameters; 
%     t=load(DUFFING_WAVEFORMS_FILE);
%     duff=t.duff_triphasic; N=2049; Fs=50; 
%
%     [sp,imp_train,est_IF]=duffing_train(duff(1,:),N,Fs);
%
%     figure(1); clf; 
%     n=(1:length(sp))./Fs; 
%     subplot(2,1,1); plot(n,sp); xlim([5 35]);
%     xlabel('time (seconds)');  ylabel('amplitude');
%     subplot(2,1,2); plot(n,est_IF); xlim([5 35]);
%     xlabel('time (seconds)');  ylabel('frequency (Hz)');
%

% John M. O' Toole, University College Cork
% Started: 05-04-2013
%-------------------------------------------------------------------------------
function [sp,imp_train,est_IF]=duffing_train(duff,N,Fs)
if(nargin<2 || isempty(N)), N=2048; end
if(nargin<3 || isempty(Fs)), Fs=50; end

DBplot=0;
DBplot_TFD=0;

STATIONARY=0;


L_d=length(duff);
duff_period=floor(L_d*0.6);
rand_limits=floor(L_d/8);

if(STATIONARY), rand_limits=0; end


%---------------------------------------------------------------------
% impulse train:
%---------------------------------------------------------------------
itrain(1)=1;  n=1;
while(itrain(n)<(N-2*(duff_period+rand_limits)))
  n=n+1;
  
  poiss_var=rand_limits*2;
  itrain(n)=itrain(n-1)+duff_period+(poissrnd(poiss_var,1,1)-poiss_var);  
end
sp=zeros(N,1);
imp_train=zeros(N,1); 
imp_train(itrain)=1;

for n=1:length(itrain)
  sp( itrain(n):(itrain(n)+L_d-1) )=duff;
end


if(length(sp)>N), sp=sp(1:N); end


Fsnew=50; 
% $$$ if(Fs/Fsnew~=1)
% $$$     sp=resample(sp,Fs/Fsnew,1); 
% $$$     N=length(sp); Fs=Fsnew;
% $$$ end


%---------------------------------------------------------------------
% estimate the IF from the impulse train
%---------------------------------------------------------------------
if(nargout>2)
  rr=diff(itrain);
  ff=1./rr;

  in=itrain(2):itrain(end);
  itrain_p=itrain(2:end);

  ff_interp=pchip(itrain_p,ff,in);    


  ff_interp_all=zeros(1,N);
  ff_interp_all(in)=ff_interp.*Fs;

  est_IF=ff_interp_all;
else
  est_IF=[];
end




%---------------------------------------------------------------------
% PLOT
%---------------------------------------------------------------------
if(DBplot)
  figure(8); clf;
  plot( (1:N)./Fs, sp);
  
  if(~isempty(est_IF))
    figure(10); clf;
    plot((1:N)./Fs, est_IF);
  end
  drawnow;
end

if(DBplot_TFD)
   dt=dtfd_sep1(sp,{11,'hamm',0,1},{1024,'bartlett'},floor(N/4),N*2); 
   figure(9); clf; 
   wide_vtfd(dt,sp,Fs,0,10);
   drawnow;
end


