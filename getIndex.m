function [pts,lPts] = getIndex(fS,iStart,iEnd,jStart,jEnd)

Nyg = fS.Nyg;

intx = iStart:iEnd;
inty = jStart:jEnd;
lIntx = (iEnd) - (iStart) + 1;
lInty = (jEnd) - (jStart) + 1;

pts = kron((intx-1)*Nyg,ones(1,lIntx)) + kron(ones(1,lInty),inty);
lPts= length(pts);

end
