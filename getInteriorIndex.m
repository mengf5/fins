function [pts,lPts] = getInteriorIndex(fS)

Nyg = fS.Nyg;
ia  = fS.ia;
ib  = fS.ib;
ja  = fS.ja;
jb  = fS.jb;

intx = ia+1:ib-1;
inty = ja+1:jb-1;
lIntx = (ib-1) - (ia+1) + 1;
lInty = (jb-1) - (ja+1) + 1;

pts = kron((intx-1)*Nyg,ones(1,lIntx)) + kron(ones(1,lInty),inty);
lPts= length(pts);

end