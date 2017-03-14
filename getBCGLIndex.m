function [pts,lPts] = getBCGLIndex(fS,axis,side,pos)
% evaluate the boundary line, ghost line indexing for given (axis,side,position)
%  position 0: boundary 1:first ghost line 2:second ghost line
Nyg = fS.Nyg;
ia  = fS.ia;
ib  = fS.ib;
ja  = fS.ja;
jb  = fS.jb;

bcxStart = (side==0)*ia + (side==1)*ib;
bcyStart = (side==0)*ja + (side==1)*jb;

bcxEnd = (axis==1)*((side==1)*ia + (side==0)*ib) +  (axis==0)*(bcxStart);
bcyEnd = (axis==0)*((side==1)*ja + (side==0)*jb) +  (axis==1)*(bcyStart);

lBcx = abs(bcxEnd - bcxStart) + 1;
lBcy = abs(bcyEnd - bcyStart) + 1;

bcx = min(bcxStart,bcxEnd):max(bcxStart,bcxEnd);
bcy = min(bcyStart,bcyEnd):max(bcyStart,bcyEnd);
% shift from boundary to ghost line wrt its pos(ition)
bcx = bcx + (axis==0)*( - (-1)^side * pos );
bcy = bcy + (axis==1)*( - (-1)^side * pos );

pts = kron((bcx-1)*Nyg,ones(1,lBcy)) + kron(ones(1,lBcx),bcy);
lPts= length(pts);


end