%-------------------------------------------------------------------------------
% thres: set to zero any value below L
%
% Syntax: x_lim=thres(x,lim)
%
% Inputs: 
%     x,lim - 
%
% Outputs: 
%     x_lim - 
%
% Example:
%     
%

% Started: 02-02-2012
%-------------------------------------------------------------------------------
function x_lim=thres(x,lim)

mSize=size(x);
x_lim=zeros(mSize);

for i=1:mSize(1)
   f=find((x(i,:))>=lim);
   x_lim(i,f)=x(i,f);
end
