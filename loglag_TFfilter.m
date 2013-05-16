%-------------------------------------------------------------------------------
% loglag_TFfilter: Filter in the log-lag domain (to remove spectral modulation)
%
% Syntax: tfd_filt=loglag_TFfilter(tfd,Fs,N,win_param)
%
% Inputs: 
%     tfd       - time-frequency distribution
%     Fs        - sampling distribution
%     win_param - parameter for Tukey window (default 0.0025)
% 
%
% Outputs: 
%     tfd_filt - filtered time-frequency distribution
%
% Example:
%      b=load('PLED_example_epoch.mat');
%      x=b.x(1:1024);
%      tf=gen_TFD_EEG(x,b.Fs,512,'sep');
%      tf_filt=loglag_TFfilter(tf,b.Fs);
%
%      figure(1); clf;
%      wide_vtfd(tf_filt,x,b.Fs,0.5,10);

% John M. O' Toole, University of Deusto
% Started: 08-05-2012
%-------------------------------------------------------------------------------
function [tfd_filt]=loglag_TFfilter(tfd,Fs,win_param)
if(nargin<2 || isempty(Fs)), Fs=1; end
if(nargin<3 || isempty(win_param)) win_param=[]; end


DB=0;
DBverbose=0;


[Ntime,M]=size(tfd);

%---------------------------------------------------------------------
% 0. parameters
%---------------------------------------------------------------------
LAG_WIN_TYPE='tukey';
LAG_WIN_LENGTH=M;
if(isempty(win_param))
    LAG_WIN_PARAM=0.0025;
else
    LAG_WIN_PARAM=win_param;
end


if(DBverbose)
    dispVars(Ntime,N,M,Fs,LAG_WIN_PARAM);
end


%---------------------------------------------------------------------
% 1. Generate filter 
%---------------------------------------------------------------------
Mh=M/2;

filt_lag=get_window(LAG_WIN_LENGTH,'tukey',LAG_WIN_PARAM);
filt_lag=shiftWin(filt_lag);

if(length(filt_lag)<M);
    filt_lag=padWin( shiftWin(filt_lag), M );
    filt_lag=shiftWin(filt_lag);
end
filt_shift=zeros(M,1);
filt_shift(1:Mh+1)=filt_lag(1:Mh+1);
filt_shift(M:-1:Mh+2)=filt_lag(2:Mh);

if(DB)
    dispVars(max(abs(imag(fft(filt_shift)))));
end

n=1:Ntime;
filt_mat=repmat(filt_shift,1,Ntime).';


%---------------------------------------------------------------------
% 2. Do the filtering in a 'cepstrum-like' domain
%---------------------------------------------------------------------
tfd_filt=ilog_DFT( ((log_IDFT(tfd.').').*filt_mat).' ).';


tfd_filt=real(tfd_filt);




function c=log_IDFT(x)
%---------------------------------------------------------------------
% Can use different log values here
%---------------------------------------------------------------------
c=ifft( log( abs(x)+eps ) );



function c=ilog_DFT(x)
%---------------------------------------------------------------------
% inverse of above
%---------------------------------------------------------------------
c=exp( fft(x) );



