%-------------------------------------------------------------------------------
% plot_results_SNR_ifests: Plot IF estimation error for the three methods
%
% Syntax: [stats_st]=plot_results_SNR_ifests(snrs)
%
% Inputs: 
%     err_st - mean-squared error (structure) for three methods
%
% Outputs: 
%     stats_st - results for t-tests
%
% Example:
%     see run_all_script_synthetic_dataset.m for examples
%

% John M. O' Toole, University College Cork
% Started: 21-05-2013
%-------------------------------------------------------------------------------
function [stats_st]=plot_results_SNR_ifests(err_st)

t_parameters;

QUARTILES=[0.25 0.5 0.75];
LOG_TRANSFORM=1;

% parameters for plotting only:
BARS_COLOR={'b','r','g'};
BARS_COLOR={[0 0 .5],'r',[0 0.5 0],[.75 .75 1]};
BARS_COLOR={[128 0 128]./256,[173 216 230]./256,[255 165 0]./256};
BARS_WIDTH=2.0;
MED_LINEWIDTH=1.2;
MED_LINE_STYLE={'-','--',':','-.'};
MARKER_SIZE=[7,7,7]-1;
MARKERS={'s','o','d'};
MARKERS_SIG={'p','*'};
MARKER_SIZE_SIG=6;
MARKER_COLOR_SIG='k';
SNR_SHIFT=0.2;



L=length(err_st);

figure(71); clf; hold all;
med_tfd=zeros(L,1); med_tfdonly=zeros(L,1); med_spike=zeros(L,1);
ip1=0; ip2=0; 

snrs=[err_st.SNR];

for n=1:L
    SNR=err_st(n).SNR;
    
    m_tfd=[]; m_tfdonly=[]; m_spike=[];
    
    m_tfd=transform_data(err_st(n).method_tfd,LOG_TRANSFORM);
    m_tfdonly=transform_data(err_st(n).method_tfdonly,LOG_TRANSFORM);    
    m_spike=transform_data(err_st(n).method_spike,LOG_TRANSFORM);        
    
    quarts_tfd=quantile(m_tfd,QUARTILES);
    quarts_tfdonly=quantile(m_tfdonly,QUARTILES);    
    quarts_spike=quantile(m_spike,QUARTILES);  
    
    med_tfd(n)=quarts_tfd(2);  
    med_tfdonly(n)=quarts_tfdonly(2);
    med_spike(n)=quarts_spike(2);    

    
    %---------------------------------------------------------------------
    % plot:
    %---------------------------------------------------------------------
    line([SNR-SNR_SHIFT SNR-SNR_SHIFT], [quarts_tfd(1) quarts_tfd(3)], ...
         'color',BARS_COLOR{1},'linewidth',BARS_WIDTH);
    line([SNR SNR], [quarts_tfdonly(1) quarts_tfdonly(3)], ...
         'color',BARS_COLOR{2},'linewidth',BARS_WIDTH);
    line([SNR+SNR_SHIFT SNR+SNR_SHIFT], [quarts_spike(1) quarts_spike(3)], ...
         'color',BARS_COLOR{3},'linewidth',BARS_WIDTH);

    
    %---------------------------------------------------------------------
    % paired t-tests:
    %---------------------------------------------------------------------
    d_tfd_tfdonly=m_tfd-m_tfdonly;
    d_tfd_spike=m_tfd-m_spike;    
    
    [p_tfdonly,h_tfdonly]=signrank(d_tfd_tfdonly,0,'tail','left');
    [p_spike,h_spike]=signrank(d_tfd_spike,0,'tail','left');    
    
    dispVars(p_tfdonly,p_spike,h_tfdonly,h_spike);
    
    if(p_tfdonly<0.0001)
        plot(SNR,1.3,'Marker',MARKERS_SIG{1}, ...
              'MarkerSize',MARKER_SIZE_SIG, ...
              'MarkerEdgeColor','k', ...
              'MarkerFaceColor',MARKER_COLOR_SIG);
    end
    if(p_spike<0.0001)
        plot(SNR,1.7,'Marker',MARKERS_SIG{2}, ...
              'MarkerSize',MARKER_SIZE_SIG, ...
              'MarkerEdgeColor','k', ...
              'MarkerFaceColor',MARKER_COLOR_SIG);
    end
    
        
end

hl=plot(snrs-SNR_SHIFT,med_tfd,snrs,med_tfdonly,snrs+SNR_SHIFT,med_spike);

for p=1:3
    set(hl(p),'linewidth',MED_LINEWIDTH, ...
              'color',BARS_COLOR{p}, ...                 
              'lineStyle',MED_LINE_STYLE{p}, ...
              'Marker',MARKERS{p}, ...
              'MarkerSize',MARKER_SIZE(p), ...
              'MarkerEdgeColor','k', ...
              'MarkerFaceColor',BARS_COLOR{p} ...                 
              );
end
hleg=legend(hl,'proposed','TFD only','spike detection');
legend boxoff;
set(hleg,'location','southwest');

xlim([min(snrs)+sign(min(snrs)) max(snrs)+sign(max(snrs))]);

xl=xlim;

ylim([-4.5 2]);
set(gca,'ytick',[-4:3]);
set(gca,'yticklabel',{'-4','-3','-2','-1','0','1','',''});  
xlabel('SNR (dB)'); ylabel('log(MSE)');




function y=transform_data(x,LOG_TRANSFORM)
%---------------------------------------------------------------------
% transform the data (using log)
%---------------------------------------------------------------------
y=x;
if(LOG_TRANSFORM), y=log10(x); end



function [] = shrinkLegend (hL, shrinkFact)
 %SHRINKLEGEND - Changes LEGEND fontsize and LEGEND axes position %
 %Syntax: shrinkLegend(hL, shrinkFact)
 %
 %Inputs:
 % hL Legend axes handle
 % shrinkFact Factor by which to shrink the legend.
 % Default is 0.8
 %
 %Example: %Make fontsize and legend axes twice bigger
 % hL=legend(......);
 % shrinklegend(hL,2);
 %
 %Authors: Jim Phillips and Denis Gilbert, 03-Dec-1999
  
 if ~exist('shrinkFact','var'), shrinkFact = 0.8; end
  
 p = get(hL, 'position');
 p(3) = p(3)*shrinkFact;
 p(4) = p(4)*shrinkFact;
 set(hL,'position', p)
 ht = findobj( get(hL,'children'), 'type', 'text');
 set(ht, 'FontSize', get(ht,'FontSize')*shrinkFact)
 set(gcf,'Resizefcn','')

% And Denis has found that for printing to postscript to work properly,
% you need:
 set(gcf,'Resizefcn','')
