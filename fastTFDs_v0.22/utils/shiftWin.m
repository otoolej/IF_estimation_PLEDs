%--------------------------------------------------------------------------------
% Shift window so that  positive indices first, then followed be negative indices
%
%--------------------------------------------------------------------------------
function w=shiftWin(w)
N=length(w);
w=circshift(w(:),ceil(N/2));
