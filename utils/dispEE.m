%% Just a simple one to disp error between 
%% matrices ..
function dispEE( a, b, msg );
 
ee = a - b;
diff_re = maxn( abs(real(ee)) );
diff_im = maxn( abs(imag(ee)) );
% $$$ a_name=deblank(argn(1,:));
% $$$ b_name=deblank(argn(2,:));
a_name=deblank(inputname(1));
b_name=deblank(inputname(2));
str=['<< ' a_name ' - ' b_name ' >>'];
str=['<< ' a_name ' - ' b_name];

str=[str ' = ' num2str(diff_re) ' + j' num2str(diff_im) ' >>'];
disp(str);
return;



disp(['<< ' a_name ' - ' b_name ' >>']);
if( nargin > 2 )
  dispVars(msg, diff_re, diff_im);
else
  dispVars(diff_re, diff_im);
end

disp('<<------------->>')

