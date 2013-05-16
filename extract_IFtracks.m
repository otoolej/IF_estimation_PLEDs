%-------------------------------------------------------------------------------
% extract_IFtracks: Estimate IF componments in time-frequency distribution using the
% method from [1]. Also, include a rule-based approach for selecting only 1 IF (which
% should be the fundamental for a haromonic signal)
%
% Syntax: if_tracks=extract_IFtracks(tf,N,Fs)
%
% Inputs: 
%     tf - time-frequency distribution
%     N  - length of time-domain signal
%     Fs - sampling frequency
%
% Outputs: 
%     if_tracks - 
%
% Example:
%     
%
% [1] McAulay, R. J., & Quatieri, T. F. (1986). Speech analysis/synthesis based on a
% sinusoidal representation. IEEE Transactions on Acoustics, Speech, and Signal
% Processing, 34(4), 744–754. 

% John M. O' Toole, University of Deusto
% Started: 20-02-2012
%-------------------------------------------------------------------------------
function [if_law,if_law_samples,f_scale,t_scale,if_tracks]= ...
        extract_IFtracks(tfd,N,Fs)
if_law=[]; if_law_samples=[]; if_tracks=[];


DBplot=0;
DB=0;
DBverbose=0;

%---------------------------------------------------------------------
% 0) set parameters:
%---------------------------------------------------------------------
DELTA_SEARCH_FREQ=0.92;  % in Hz/s 
MIN_IF_LENGTH=6;      % in seconds

% limit the search in the TFD to this area:
LIMIT_TRACK_SEARCH=1;
lower_limit=0.1; upper_limit=5;

% if allowing to return a second IF:
ADD_SECOND_IF=1;
OVERLAP_ALLOWED=0.1;
IF_RANGE_WEIGTH=0.2;


[Ntime,Nfreq]=size(tfd);

dur=N/Fs;
t_scale=(dur/Ntime);
f_scale=(1/Nfreq)*(Fs/2);


nlower=ceil(lower_limit/f_scale);
nupper=ceil(upper_limit/f_scale);        

if(LIMIT_TRACK_SEARCH)
    tfd=tfd(:,nlower:nupper);    
else
    lower_limit=[];
end
    

delta_freq_samples=floor( (DELTA_SEARCH_FREQ/f_scale)*t_scale );
min_component_length=floor( MIN_IF_LENGTH/t_scale );

%---------------------------------------------------------------------
% 1) estimate IF tracks (using method in [1])
%---------------------------------------------------------------------
if_tracks=tracks_MCQmethod(tfd,Fs,delta_freq_samples,min_component_length);

  
if(isempty(if_tracks))
    return;
end
    
%---------------------------------------------------------------------
% 2) Rule-set for extracting 'main' IF
%---------------------------------------------------------------------

% a. only interested in tracks with a median value within this range:
med_lower_freq=0.2; med_upper_freq=upper_limit;        

range_tracks=[floor(med_lower_freq/f_scale) ...
              ceil(med_upper_freq/f_scale)];

Ntracks=length(if_tracks);

if_energy=[]; if_med=[]; if_var=[];

p=1;
for i=1:Ntracks
    if(size(if_tracks{i},1)>0)
        ifest=if_tracks{i}(:,2)+nlower; 
        if_median=median(ifest);

        % if between the range we are interested in
        if( if_median>=range_tracks(1) & if_median<=range_tracks(2) )
            subset_tracks{p}(:,1)=if_tracks{i}(:,1);
            subset_tracks{p}(:,2)=ifest;
            
            n=if_tracks{i}(:,1);
            k=if_tracks{i}(:,2);

            
            % features for each IF:
            if_energy(p)=0;
            for q=1:length(n)
                if_energy(p)=if_energy(p)+sum( tfd(n(q),k(q)) );
            end
            if_med(p)=if_median*f_scale;
            if_var(p)=std( ifest*f_scale );
            if_length(p)=length(n).*t_scale;
            
            p=p+1;
        end
    end
end


if(DB)
    dispVars(if_energy);
    dispVars(if_med);
    dispVars(if_var);        
end


% b. select IF with maxium energy:

[d,imax_energy]=max(if_energy);
if(isempty(if_energy) | d<0)
    if_law=[]; if_law_samples=[];
else
    if_law=zeros(length(subset_tracks{imax_energy}(:,1)),2);    
    if_law(:,1)=subset_tracks{imax_energy}(:,1).*t_scale;
    if_law(:,2)=subset_tracks{imax_energy}(:,2).*f_scale;    
    
    if_law_samples=zeros(Ntime,1);
    if_law_samples(subset_tracks{imax_energy}(:,1))= ...
        subset_tracks{imax_energy}(:,2).*f_scale;
    
    % add a rule to include a second IF:
    if(ADD_SECOND_IF)
        if_range=[if_med(imax_energy)-IF_RANGE_WEIGTH*if_med(imax_energy),
                  if_med(imax_energy)+IF_RANGE_WEIGTH*if_med(imax_energy)];
        
        % 1) find all within range (exclude one already found)
        iposs=find(if_med>if_range(1) & if_med<if_range(2));
        iposs( (iposs==imax_energy) )=[];
        
        % 2) if the overlap in time is within a limit combine two IF laws
        Lif_overlap=floor(OVERLAP_ALLOWED*length(if_law)/t_scale);
        ii=1; greater_energy=1; if_new=[];
        for i=1:length(iposs)
            itest=iposs(i); 
            if_cand=zeros(length(subset_tracks{itest}(:,1)),2);    
            if_cand(:,1)=subset_tracks{itest}(:,1).*t_scale;
            if_cand(:,2)=subset_tracks{itest}(:,2).*f_scale;    
            
            [c,ic,il]=intersect(if_cand(:,1),if_law(:,1));
% $$$                 dispVars(ii,length(c),Lif_overlap);
            if(length(c)<Lif_overlap)
                if(ii>1)
                    if(if_energy(ii)<if_energy(ii-1))
                        greater_energy=0;
                    end
                end
                if(greater_energy==1 & if_law(1,1)<=if_cand(1,1))
                    new_time=cat(1,if_law(:,1), if_cand(length(c)+1:end,1));
                    new_freq=cat(1,if_law(:,2), if_cand(length(c)+1:end,2));

                    if_new=zeros(length(new_time),2);
                    if_new(:,1)=new_time;
                    if_new(:,2)=new_freq;                        
                    
                elseif(greater_energy==1 & if_law(1,1)>if_cand(1,1))
                    new_time=cat(1,if_cand(1:end-length(c)-1,1),if_law(:,1));
                    new_freq=cat(1,if_cand(1:end-length(c)-1,2),if_law(:,2)); 
                    
                    if_new=zeros(length(new_time),2);
                    if_new(:,1)=new_time;
                    if_new(:,2)=new_freq;                        
                end
                ii=ii+1;
            end
            greater_energy=1;
        end
        if(~isempty(if_new))
             if(DBverbose) dispVars('---- new IF law ----'); end
            if_law=if_new;
        end
    end
end

% paste if_law into to if_law_samples aswell:
if_law_samples=zeros(Ntime,1);
if(~isempty(if_law))
    tt=round(if_law(:,1)./t_scale);
    if_law_samples(tt)=if_law(:,2)./f_scale;
end



if(DBplot)
    figure(65); clf;
    if(~isempty(if_law))
        plot(if_law(:,1),if_law(:,2),'+'); hold on;
        plot((1:Ntime).*t_scale,if_law_samples.*f_scale,'r');
    end
    axis([0 dur lower_limit upper_limit]);
end


