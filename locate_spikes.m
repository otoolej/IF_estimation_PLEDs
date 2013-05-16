%-------------------------------------------------------------------------------
% locate_enhance: Methods (from [1]-[3]) to find spikes;
%
% Syntax: [sp_enhanc,sp_train]=locate_spikes(x,Fs,METHOD,HP_FILTER)
%
% Inputs: 
%     x         - time-domain signal
%     Fs        - sampling frequency
%     METHOD    - either: 'NEO' from [1], 'SARKELA' from [2], or 'DEBURCH' from [3]
%     HP_FILTER - high-pass filter (3 Hz), 1=yes, 0=no
%
% Outputs: 
%     sp_enhanc - enhanced spikes without thresholding
%     sp_train  - enhanced spikes with thresholding
%
% Example:
%      b=load('PLED_example_epoch.mat');
%      [se,st]=locate_spikes(b.x,b.Fs,'DEBURCH',0);
%      L=length(b.x);
%      figure(1); clf; plot(1:L,b.x,1:L,st);
%      figure(2); clf; plot(1:L,b.x,1:L,se);
%      
%     
% [1] Kaiser, J. On a simple algorithm to calculate the ‘energy’ of a signal.  In:
% Int. Conf. Acoustics, Speech, and Signal Process., ICASSP-90. IEEE; 1990. p. 381–384.
%
% [2] Särkelä M, Mustola S, Seppänen T, Koskinen M, Lepola P, Suominen K, Juvonen T,
% Tolvanen-Laakso H, Jäntti V. Automatic analysis and monitoring of burst suppression in
% anesthesia. J Clin Monit Comput 2002; 17:125–34
% 
% [3] Deburchgraeve, W., Cherian, P. J., De Vos, M., Swarte, R. M., Blok, J. H., Visser,
% G. H., Govaert, P., et al. (2008). Automated neonatal seizure detection mimicking a
% human observer reading EEG. Clinical Neurophysiology, 119(11), 2447–2454. 
% 


% John M. O' Toole, University College Cork
% Started: 09-05-2013
%-------------------------------------------------------------------------------
function [sp_enhanc,sp_train]=locate_spikes(x,Fs,METHOD,HP_FILTER)
if(nargin<2 || isempty(Fs)), Fs=1; end
if(nargin<3 || isempty(METHOD)), METHOD='NEO'; end
if(nargin<4 || isempty(HP_FILTER)), HP_FILTER=1; end

x=real(x);
N=length(x);
DB=0;

%---------------------------------------------------------------------
% 1) Filter:
%---------------------------------------------------------------------
if(HP_FILTER)
    HP_fc=3;
    L_hpf=241;

    x=filter_zerophase(x,Fs,[],HP_fc,L_hpf);
end



%---------------------------------------------------------------------
% 2) emphasize spikes
%---------------------------------------------------------------------
switch METHOD
    
  %-----------------------------------------------------------
  % NLEO, from [3]
  %-----------------------------------------------------------
 case 'SARKEL'
  n=3:(N-1); 
  l=0; p=3; q=1; s=2;
  neo(n+1)=x(n-l+1).*x(n-p+1) - x(n-q+1).*x(n-s+1);

  % Get peaks: 
  Nneo=length(neo);
  sp_train=fpeaks(neo,5*median(neo));    

  sp_enhanc=neo;
  
  
  %-----------------------------------------------------------
  % NEO (nonlinear energy operator), from [1]
  %-----------------------------------------------------------
 case 'NEO'
  n=1:(N-2); 
  neo(n+1)=x(n+1).^2 - x(n+2).*x(n);
  
  % Get peaks: 
  Nneo=length(neo);
  sp_train=fpeaks(neo,5*median(neo));    

  sp_enhanc=neo;

  
  %-----------------------------------------------------------
  % Smoothed NLEO, from [2]
  %-----------------------------------------------------------
 case 'DEBURCH'
  n=3:(N-1); 
  l=1; p=2; q=0; s=3;
  neo(n+1)=x(n-l+1).*x(n-p+1) - x(n-q+1).*x(n-s+1);

  % Get peaks: 
  Nneo=length(neo);
  
  % Filter NEO: (moving average with 120ms window)
  Lw=floor( 0.120*Fs );
  neo=filter( ones(Lw,1)./Lw,1,neo); 
  
% $$$   thres=0.6*(std(neo)+percentile(neo,3));
% matlab function is 'prctile'
  thres=0.6*(std(neo)+prctile(neo,75));  

  sp_train=fpeaks(neo,thres);    

  sp_enhanc=neo;  


 otherwise
  error('Which type of method?')
end


if(nargout>2)
  Nsp=min( length(sp_train), length(sp_enhanc) );
end






function [pe,xe]=fpeaks(x,thres)
%---------------------------------------------------------------------
% Find the maximum points (peaks) above certain threshold
%
% (should vectorize this function to improve speed)
%---------------------------------------------------------------------
N=length(x);
xe=ones(N,1)*thres; pe=ones(N,1)*thres;
b=find(x>=thres);
xe(b)=x(b);
i_thres=find(xe==thres);
DB=0;


Nc=length(i_thres);

segment=xe(1:i_thres(1));
i_mx=find(segment==max(segment));
pe(i_mx)=segment(i_mx);

for n=1:Nc-1
  if(DB) dispVars('< --- >',xe(n),i_thres(n));  end

  if( i_thres(n+1)-i_thres(n)>1 )
    i_seg=[i_thres(n):i_thres(n+1)];
    segment=xe(i_seg);
    if(DB) dispVars(segment,i_seg); end
    i_mx=find(segment==max(segment));
    i_pe=i_mx+i_thres(n)-1;
    pe(i_pe)=segment(i_mx);
    if(DB) 
      dispVars(i_thres(n+1),i_mx,i_pe,max(segment),pe(i_pe)); 
      pause;
    end    
  end
end
  
