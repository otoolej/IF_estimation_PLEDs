function y=zero_pad(x,Npad)
N=length(x);

if(size(x,1)<2) x=x.'; end
    

if(Npad<N)
    error('Npad is less than N.');
elseif(Npad==N)
    y=x;
else
    y=zeros(1,Npad);
    
    y(1:N)=x;
    y(N+1:Npad)=0;    
end
    
