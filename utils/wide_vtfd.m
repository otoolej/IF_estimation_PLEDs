%-------------------------------------------------------------------------------
% Syntax: wide_vtfd(tfd,s1,FS,f_bottom,f_top,PRINT,fname,PLOT_SIZE,st_fonts)
%
% Plot TFD and time/frequency plots
%
% TFD can be decimated; will then get time info from time-domain
% signal s1
%
% Inputs: 
%     tfd,s1,FS,f_bottom,f_top,PRINT,fname,PLOT_SIZE,st_fonts - 
%
% Outputs: 
%      - 
%
% Example:
%     
%  STARTED: 08-09-2010

%   Copyright (c) 2010, John M. O' Toole, The University of Queensland
%   All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following
%  conditions are met:
%      * Redistributions of source code must retain the above
%        copyright notice, this list of conditions and the following
%        disclaimer.
%      * Redistributions in binary form must reproduce the above
%        copyright notice, this list of conditions and the following
%        disclaimer in the documentation and/or other materials
%        provided with the distribution.
%      * Neither the name of the <organization> nor the names of its
%        contributors may be used to endorse or promote products
%        derived from this software without specific prior written
%        permission.
%  
%  THIS SOFTWARE IS PROVIDED BY JOHN M. O' TOOLE ''AS IS'' AND ANY
%  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
%  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL JOHN M. O' TOOLE BE
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
%  OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
%  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
%  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
%  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
%  DAMAGE.
%-------------------------------------------------------------------------------
function wide_vtfd(tfd,s1,FS,f_bottom,f_top,PRINT,fname,PLOT_SIZE,st_fonts)
if(nargin<8 || isempty(PLOT_SIZE)) PLOT_SIZE='long'; end
if(nargin<9 || isempty(st_fonts)) st_fonts=[]; end





if(nargin<3) FS=1; end
if(nargin<2 || isempty(s1)) 
  TFD_ONLY=1; 
  Ntime=size(tfd,1);  
else 
  TFD_ONLY=0; 
  Ntime=length(s1);
end
if(nargin<6 || isempty(PRINT)) PRINT=0; end
if(nargin<7 || isempty(fname)) fname='test.eps'; end

[N,M]=size(tfd);

% If want to plot only a segment of TFD (in the frequency dir.)
SET_F_RANGE=0;
if(nargin<4 || isempty(f_bottom))
    f_bottom=0;
else
    SET_F_RANGE=1;
end
if(nargin<5 || isempty(f_top))
    f_top=(FS/2);
else
    SET_F_RANGE=1;
end




if(PRINT) 
    COLOR_BAR=0;

    if(isempty(st_fonts))
        FONT_TYPE='Times-Roman';
    
        switch PLOT_SIZE
          case 'long'
            FONT_SIZE=10;
          case 'medium'
            FONT_SIZE=16;    
          case 'short'
            FONT_SIZE=22;
          otherwise
            FONT_SIZE=10;
        end
    else
        FONT_TYPE=st_fonts{1};
        FONT_SIZE=st_fonts{2};
    end
    
else
    COLOR_BAR=1; 
end


if(iscell(s1))
  Na=length(s1);
  if(Na<2)
    error('Need two signals if cell.');
  end
  x=s1{1}; y=s1{2};
  x=x(:); y=y(:);
  Ntime=length(x);
  if(length(x)~=length(y))
    error('Time-domain signals must be the same size.');
  end
  
  S1=fft(x);
  MULTI_TIME_SIG=1;
else
  MULTI_TIME_SIG=0;
  S1=fft(s1);
end

clf;


%---------------------------------------------------------------------
% Frame sizes:
%---------------------------------------------------------------------
if(COLOR_BAR) 
  X_CBAR_START=0.03; 
  X_CBAR_WIDTH=0.025;
  X_COLOR_BAR_GAP=0.03;
else
  X_COLOR_BAR_GAP=0; 
end

if(PRINT)
    SHIFT_TF_IN_Y=0.05;
    SHIFT_TF_IN_Y=0.15;    
else
    SHIFT_TF_IN_Y=0;
end

% HACK, but may want to p this. Loo better on Alvaro's PC
SHIFT_TF_IN_Y=0.1;   



TIME_FREQ_PLOTS_WIDTH=0.1;
TIME_FREQ_PLOTS_GAP=0.07;
% $$$ TIME_FREQ_PLOTS_GAP=0.1;

X_AXIS_WIDTH=0.83-X_COLOR_BAR_GAP;
Y_AXIS_HEIGHT=0.75-SHIFT_TF_IN_Y;    

X_AXIS_START=0.01+X_COLOR_BAR_GAP;
Y_AXIS_START=0.23+SHIFT_TF_IN_Y;



%---------------------------------------------------------------------
% X/Y tick labels
%---------------------------------------------------------------------
ntime=1:Ntime; ntime=ntime./FS;
n=linspace(ntime(1),ntime(end),N);
Mh=ceil(M/2);
Ts=1/FS;

% this for spectrum:
kbottom=fix(f_bottom*Ntime*Ts);  ktop=fix(f_top*Ntime*Ts);
k=(kbottom+1:ktop);
N_k=length(kbottom:ktop);
% this for TFD:  (assming here that frequency domain is sampled with k/2NT)
if(SET_F_RANGE)
    kbottom=fix(f_bottom*2*M*Ts);  ktop=fix(f_top*2*M*Ts);
    k_tfd=(kbottom+1:ktop);
    f_tfd=k_tfd./(2*M*Ts);
else
    k_tfd=1:M;
    f_tfd=linspace(0,FS/2,M);
end




if(TFD_ONLY)
  imagesc(k,n,tfd); axis xy;
  
else
  %---------------------------------------------------------------------
  % 1. freq plot
  %---------------------------------------------------------------------
  h_freq=subplot(2,2,2);
  set(h_freq,'position',[(X_AXIS_WIDTH+TIME_FREQ_PLOTS_GAP+X_COLOR_BAR_GAP) Y_AXIS_START ...
                      abs(TIME_FREQ_PLOTS_WIDTH-X_AXIS_START+0.5*TIME_FREQ_PLOTS_GAP) Y_AXIS_HEIGHT]);
  p_freq=plot(abs(S1(k)).^2,k);  
  axis('tight');
  set(h_freq,'xticklabel',[]); set(h_freq,'yticklabel',[]);
  set(h_freq,'xtick',[]);   
  if(PRINT) 
      set(p_freq,'Linewidth',2);
  else
% $$$       grid('on'); 
  end


  %---------------------------------------------------------------------
  % 2. time plot
  %---------------------------------------------------------------------
  h_time=subplot(2,2,3);
  set(h_time,'position',[X_AXIS_START ...
                      (Y_AXIS_START-2*TIME_FREQ_PLOTS_GAP-TIME_FREQ_PLOTS_WIDTH) ...
                      X_AXIS_WIDTH TIME_FREQ_PLOTS_WIDTH]);
  if(MULTI_TIME_SIG)
    p_time=plot(ntime,x,'b',ntime,y.*max(x),'r');  
  else
    p_time=plot(ntime,s1);  
  end
  axis('tight'); 
  set(h_time,'xticklabel',[]); set(h_time,'yticklabel',[]);
  set(h_time,'ytick',[]);
  if(PRINT)
      set(p_time,'Linewidth',2);
  else
      grid('on');
  end
  


  %---------------------------------------------------------------------
  % 2. time-freq plot
  %---------------------------------------------------------------------
  h_image=subplot(2,2,1);
  set(h_image,'position',[X_AXIS_START Y_AXIS_START X_AXIS_WIDTH ...
                      Y_AXIS_HEIGHT]);
  imagesc(n,f_tfd,tfd(:,k_tfd).'); axis xy;
  set(h_image,'YAxisLocation','right');
  
  if(COLOR_BAR) 
    h_colorb=colorbar('west');
    set(h_colorb,'Position',[X_CBAR_START Y_AXIS_START X_CBAR_WIDTH ...
                        Y_AXIS_HEIGHT]);
    set(h_colorb,'YAxisLocation','left');  
  end
  
  if(PRINT)
      xlabel('TIME (seconds)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);
      ylabel('FREQUENCY (Hz)','FontName',FONT_TYPE,'FontSize',FONT_SIZE);      
      set(gca,'FontName',FONT_TYPE,'FontSize',FONT_SIZE);      
  end

end


if(PRINT)
    % turn off ticks ?
    set(h_time,'xtick',[]); 
    set(h_freq,'ytick',[]);    
    
    switch PLOT_SIZE
      case 'long'
        SIZE_PLOT='-S650,150';
      case 'medium'
        SIZE_PLOT='-S650,300';
      case 'short'
        SIZE_PLOT='-S600,550';
      otherwise
        SIZE_PLOT='-S650,150';
    end
    
    % invert to gray colormap:
% $$$     colormap( flipud(bone) );
    
    
    colorbar('off');
    ffonts=['-F', FONT_TYPE, ':', num2str(FONT_SIZE) ];
    print(fname,'-depsc2',SIZE_PLOT,ffonts);
    system( ['gv ', fname, ' &'] );
    
    system( ['convert_eps2pdf.bash ', fname, ' &'] );

end

