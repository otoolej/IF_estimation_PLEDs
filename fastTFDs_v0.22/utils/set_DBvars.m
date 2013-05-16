%-------------------------------------------------------------------------------
% FUNCTION to set global DEBUGGING variables, defined in global_vars.m
%

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
function set_DBvars(TYPE,value)
global FAST_DTFD_DB;
global FAST_DTFD_DBvars;
global FAST_DTFD_DBwarn;
global FAST_DTFD_DBplot;
global FAST_DTFD_DBtfd;
global FAST_DTFD_DBtest;

if(nargin<2) error('Need two in-args here.'); end
TYPE=lower(TYPE);


switch TYPE
 case 'all'
   FAST_DTFD_DB=value;
   FAST_DTFD_DBwarn=value;
   FAST_DTFD_DBvars=value;
   FAST_DTFD_DBplot=value;
   FAST_DTFD_DBtfd=value;
   FAST_DTFD_DBtest=value;
   FAST_DTFD_DBprofile=value;   

 case 'db'
   FAST_DTFD_DB=value;

 case 'dbvars'
   FAST_DTFD_DBvars=value;

 case 'dbplot'
   FAST_DTFD_DBplot=value;
   
 case 'dbtfd'
   FAST_DTFD_DBtfd=value;

 case 'dbwarn'
   FAST_DTFD_DBwarn=value;
   
 case 'dbtest'
   FAST_DTFD_DBtest=value;
   
 case 'dbprofile'
  FAST_DTFD_DBprofile=value;
   
   
 otherwise
  error('Incorrect DB value to set.');
end




