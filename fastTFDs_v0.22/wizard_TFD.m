%-------------------------------------------------------------------------------
% Wizard to help generate a time-frequency distribution (TFD)
% 
% USE: tf=wizard_TFD(x)
%
% INPUT:
%      x = time-domain signal (optional)
%
% OUTPUT:
%      tf = TFD 
%

%  Copyright (c) 2010, John M. O' Toole, The University of Queensland
%  All rights reserved.
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
%      * Neither the name of the The University of Queensland nor the 
%        names of its contributors may be used to endorse or promote 
%        products derived from this software without specific prior 
%        written permission.
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
function tf=wizard_TFD(x)
if(nargin<1 || isempty(x)) x=[]; end

DB=0;

win_dopp=[]; win_lag=[];
time_dec=[]; freq_dec=[];
Nfreq=[]; Ntime=[]; tf=[];


disp_sep('MESSAGE: Welcome');
disp('This function will guide you to compute');
disp('a time-frequency distribution (TFD)');


%---------------------------------------------------------------------
% 1. Time-domain signal
%---------------------------------------------------------------------
if(isempty(x))
    disp_sep('INPUT: Time-domain Signal');
     disp(['Enter time-domain signal as an already']);
    disp(['assigned variable or a MATLAB expression.']);
    x=input('default value gen_LFM(128,0.1,0.3): ');
    if(isempty(x))
        x=gen_LFM(128,0.1,0.4);
    end
 end
if(DB) dispVars(size(x)); end
N=length(x);
    
    
%---------------------------------------------------------------------
% 2. Type of TFD
%---------------------------------------------------------------------
disp_sep('INPUT: TFD Type');
kern_type=list_menu('kernel types:', ...
                   'Smoothed-Pseudo Wigner-Ville (separable kernel)', ...           
                   'Choi-Williams (nonseparable kernel)', ...
                   'Pseudo Wigner-Ville (Doppler-independent kernel)', ...
                   'Smoothed Wigner-Ville (lag-independent kernel)' );
% $$$                    'Wigner-Ville');
if(DB) dispVars(kern_type); end



%---------------------------------------------------------------------
% 3. Kernel type
%---------------------------------------------------------------------
if(kern_type~=5)
    disp_sep('INPUT: Kernel Parameters');
    
    switch kern_type
      % nonsep kernel:
      case 1
        kern_types='cw';
        disp('Enter kernel parameter for Choi-Williams distribution.');
        kern_param=pick_integer_number(10);
    
        % Sep kernel:
      case 0
        kern_types='sep';
        disp_sep('INPUT: Doppler Window');
        win_dopp=window_params(N,1);
        disp_sep('INPUT: Lag Window');
        win_lag=window_params(N,0);
        kern_param={win_dopp,win_lag};
        
        % DI kernel:
     case 2
        kern_types='pwvd';
        disp_sep('INPUT: Lag Window');
        win_lag=window_params(N,0);
        kern_param=win_lag;
      
        % LI kernel:
     case 3
        kern_types='swvd';
        disp_sep('INPUT: Doppler Window');
        win_dopp=window_params(N,1);
        kern_param=win_dopp;
      
     otherwise
      error('which kernel type?');
    end
    
    % special case for the WVD distribution:
elseif(kern_type==5)
  kern_type='wvd';
  kern_params={};
end
if(DB) dispVars(kern_types); end



%---------------------------------------------------------------------
% 4. Nfreq and Ntime frequency and time-domain sample rate
%---------------------------------------------------------------------
if(~isempty(win_lag))
  Nfreq_min=ceil( (win_lag{1}+1)/2 );
  Nfreq_default=2^(nextpow2(Nfreq_min));
  
  disp_sep('INPUT: length of TFD in frequency direction');
  Nfreq=pick_integer_number(Nfreq_default,Nfreq_min);
end
if(~isempty(win_dopp))
  Ntime_min=2*win_dopp{1};
  Ntime_default=2^(nextpow2(Ntime_min));
  
  disp_sep('INPUT: length of TFD in time direction');  
  Ntime=pick_integer_number(Ntime_default,Ntime_min);
end


%---------------------------------------------------------------------
% 5. Decimate?
%---------------------------------------------------------------------
disp_sep('INPUT: decimated TFD'); 
yn=input('Produce a decimated TFD? Input yes/no [no]: ','s');
if(isempty(yn)) yn='no'; end % redundant code for octave
if(strcmp(yn(1),'y'))

% $$$   disp_sep('INPUT: length of TFD in time direction');  
% $$$   Ntime=pick_integer_number(Ntime_default,Ntime_min)

  % DI kernel:
  if(kern_type==2)
    
    time_dec_ok=0;
    while(time_dec_ok==0)
      disp(['Enter either a scaler value or a vector [n1,n2,...,nU]' , ...
            'where U<2N  and 0<=ni<=2N (N is signal length)']);
      time_dec=input('Input time decimatation factor [2]: ');

      time_dec=fix(time_dec);
      if(length(time_dec)>1)
        iout_of_range=find(time_dec>2*N || time_dec<1);
        if(isempty(iout_of_range))
          time_dec_ok=1;
        end
      else
        % not checking if scalar value is ok
        time_dec_ok=1;
      end
    end

    disp('Frequency decimation scalar:');
    freq_dec=pick_integer_number(2);
  end

  % LI kernel:
  if(kern_type==3)
    
    freq_dec_ok=0;
    while(freq_dec_ok==0)
      disp(['Enter either a scaler value or a vector [n1,n2,...,nV]' , ...
            'where V<N  and 0<=ni<=N (N is signal length)']);
      freq_dec=input('Input frequency decimatation factor [2]: ');

      freq_dec=fix(freq_dec);
      if(length(freq_dec)>1)
        iout_of_range=find(freq_dec>N || freq_dec<1);
        if(isempty(iout_of_range))
          freq_dec_ok=1;
        end
      else
        % not checking if scalar value is ok
        freq_dec_ok=1;
      end
    end

    disp('Time decimation scalar:');
    time_dec=pick_integer_number(2);
  end
  
  % nonsep. or sep. kernel:
  if(kern_type==1 || kern_type==0)
      disp('Time decimation scalar:');
      time_dec=pick_integer_number(2);
      disp(' ');
      disp('Frequency decimation scalar:');
      freq_dec=pick_integer_number(2);
  end
  
end
if(DB) dispVars(time_dec,freq_dec); end
    

%---------------------------------------------------------------------
% 6. Assemble function
%---------------------------------------------------------------------

% if sep. kernel:
if(kern_type==0)
    win_lag_str=win_cell_to_string(win_lag);
    win_dopp_str=win_cell_to_string(win_dopp);
    kern_param_str=['{' win_dopp_str ',' win_lag_str '}'];
    
    % or nonsep. kernel:
elseif(kern_type==1)
    kern_param_str=['{', num2str(kern_param), '}'];
    
    % or DI or LI kernel:
else
    kern_param_str=win_cell_to_string( kern_param );
end

if(isempty(Nfreq))   Nfreq=N; end
if(isempty(Ntime))   Ntime=2*N; end
Nfreq_str=num2str(Nfreq);
Ntime_str=num2str(Ntime);

if(isempty(time_dec))
    time_dec_str='no';
else
    if(length(freq_dec)>1)
        time_dec_str=['[' num2str(time_dec) ']'];
    else
        time_dec_str=num2str(time_dec);
    end
end
if(isempty(freq_dec))
    freq_dec_str='no';
else
    if(length(freq_dec)>1)
        freq_dec_str=['[' num2str(freq_dec) ']'];
    else
        freq_dec_str=num2str(freq_dec);
    end
    
end

    
clc;
disp_sep('MESSAGE:  Compute TFD'); 

disp(['Kernel type                           : ', kern_types]);
disp(['Kernel parameter                      : ', kern_param_str]);
disp(['Signal length                         : ', num2str(length(x))]);
disp(['Nfreq (dimension in freq. direction ) : ' Nfreq_str]);
disp(['Ntime (dimension in time direction )  : ' Ntime_str]);
disp(['Time decimation                       : ' time_dec_str]);
disp(['Frequency decimation                  : ' freq_dec_str]);


an=input('Proceed to compute TFD? yes/no [yes]: ','s');
if(isempty(an)) an='yes'; end  % redundant code for octave

if(an(1)=='y')
    
    % if FULL (or exact) TFD:
    if( isempty(time_dec) && isempty(freq_dec) )
        switch kern_type
            % nonsep kernel:
          case 1
            tf=dtfd_nonsep(x,kern_types,{kern_param});
            % Sep kernel:
          case 0
            tf=dtfd_sep1(x,win_dopp,win_lag,Ntime,Nfreq);
            % DI kernel:
          case 2
            tf=dtfd_DI(x,kern_param,Nfreq);
            % LI kernel:
          case 3
            tf=dtfd_LI(x,kern_param,Ntime);
          otherwise
            error('which kernel type?');
        end

        
        % otherwise decimated  TFD:
    else
        switch kern_type
            % nonsep kernel:
          case 1
            tf=dec_dtfd_nonsep(x,kern_types,{kern_param},time_dec,freq_dec);
            % Sep kernel:
          case 0
            tf=dec_dtfd_sep(x,win_dopp,win_lag,Ntime,Nfreq,time_dec,freq_dec);
            % DI kernel:
          case 2
            tf=dec_dtfd_DI(x,kern_param,Nfreq,time_dec,freq_dec);
            % LI kernel:
          case 3
            tf=dec_dtfd_LI(x,kern_param,Ntime,time_dec,freq_dec);
          otherwise
            error('which kernel type?');
        end

    end    
    
    disp_sep('MESSAGE: Plot TFD');
    an=input('Plot TFD? yes/no [yes]: ','s');
    if(isempty(an)) an='yes'; end  % redundant code for octave
    if(an(1)=='y')
        clf; vtfd(tf,x);
    end
    
end





function disp_sep(msg)
disp(' ');
disp(['--------------------  ' msg '  --------------------' ]);
disp(' ');



function in=list_menu(heading,varargin)
%---------------------------------------------------------------------
% simple text menu 
%---------------------------------------------------------------------
Nitems=length(varargin);
in_ok=0;

% 1. draw menu:
while(in_ok==0)
    disp(sprintf(heading));
% $$$     disp(' ');

    item_nums=0:Nitems-1;
    for n=1:Nitems
        disp([' ' num2str(item_nums(n)) ') ' sprintf(varargin{n})]);
    end
    disp(' ');

    % 2. get input
    in=input('Select a menu number [default 0]: ');
    if(isempty(in)) in=0; end


    % 3. check if value is ok?
    inum=find(item_nums==in);
    if(~isempty(inum))
        in_ok=1;
    end
end



function win_cell=window_params(N,doppler_win)
%-------------------------------------------------------------------------------
% get user to input parameters for window function
%-------------------------------------------------------------------------------
if(doppler_win)
    Ldefault=make_odd(fix(N/5));
    Lmax=make_odd(N-2);
    domain={'Doppler','time'};
else
    Ldefault=make_odd(fix(N/2));
    Lmax=2*N-1;
    domain={'lag','frequency'};
end


% 1. which domain?
idomain=list_menu('window defined in:', ...
                  [domain{1} '-domain'], ...
                  [domain{2} '-domain']);
if(doppler_win)
    win_dft=abs(idomain-1);
else
    win_dft=idomain;
end


% 2. get window length
L_ok=0;
while(L_ok==0)
    L=input(['Enter window length (integer value between 3 and ' ...
             num2str(Lmax) ')  [' num2str(Ldefault) ']:']);
    if(isempty(L)) L=Ldefault; end
    L=fix(L);
    if(L>=3 && L<=Lmax)
        L_ok=1;
    end
end
     
% 3. get window type
windows={ 'Hamming', ...
          'Hanning', ...
          'Rectangle', ...
          'Bartlett', ...          
          'Guassian', ...
          'Tukey', ...
          'Cosh' };
iwin=list_menu('window type:', ...
                  windows{1}, windows{2}, ... 
                  windows{3}, windows{4}, ...                   
                  windows{5}, windows{6}, ...                      
                  windows{7});
win_type=windows{iwin+1};


% 4. get window parameter
win_param=[];
if( strcmp(win_type,'Guassian')==1 || ...
    strcmp(win_type,'Tukey')==1 || ...
    strcmp(win_type,'Cosh')==1 )
  param_ok=0;
  while(param_ok==0)
    win_param=input('Enter window parameter:'); 
    if(isnumeric(win_param))
      param_ok=1;
    end
  end
end


win_cell{1}=L;
win_cell{2}=win_type;
if(~isempty(win_param))
  win_cell{3}=win_param;
end
if(~isempty(win_dft))
  win_cell{4}=win_dft;
end

% 5. display output:
win_str=win_cell_to_string(win_cell);
disp(' ');
disp(['Window function is {win_length,win_type,win_param,win_dft}: ' ...
      win_str]);



function win_str=win_cell_to_string(win_cell)
%-------------------------------------------------------------------------------
% Convert cell to string
%-------------------------------------------------------------------------------
if(isempty(win_cell{3}))
  win_param_str='[]';
else
  win_param_str=num2str(win_cell{3});
end

switch length(win_cell)
 case 2
  win_str=['{' num2str(win_cell{1}) ',''' win_cell{2} '''}'];
 case 3
  win_str=['{' num2str(win_cell{1}) ',''' win_cell{2} ''',' ...
           win_param_str '}'];
 case 4
  win_str=['{' num2str(win_cell{1}) ',''' win_cell{2} ''',' ...
             win_param_str ',' num2str(win_cell{4}) '}'];
end

  
function b=pick_integer_number(b_default,b_min)
%-------------------------------------------------------------------------------
% Input an integer number >b_min
%-------------------------------------------------------------------------------
if(nargin<2) b_min=[]; end
b_ok=0;

if(~isempty(b_min))
    while(b_ok==0)
        b=input(['Input integer number (min value: ' num2str(b_min) ') [' num2str(b_default) ...
                 '] :']);
  
        if(isempty(b)) b=b_default; end
  
        if(b>=b_min) b_ok=1; end
    end
else
    while(b_ok==0)
        b=input(['Input integer number [' num2str(b_default) ']: ']);
  
        if(isempty(b))
            b=b_default;
        end
        b_ok=1;
    end
end



function y=make_odd(x)
%---------------------------------------------------------------------
% return odd integer
%---------------------------------------------------------------------
if(~rem(x,2))
    y=x+1;
else
    y=x;
end
