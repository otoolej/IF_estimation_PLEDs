function x=get_odd_int(x)
x=fix(x);
if(mod(x,2))
    x=x-1;
end
