function [pts,lPts] = getIndex(fS,iStart,iEnd,jStart,jEnd)

Nyg = fS.Nyg;

intx = iStart+1:iEnd-1;
inty = jStart+1:jEnd-1;
lIntx = (ib-1) - (ia+1) + 1;
lInty = (jb-1) - (ja+1) + 1;

pts = kron((intx-1)*Nyg,ones(1,lIntx)) + kron(ones(1,lInty),inty);
lPts= length(pts);

end
