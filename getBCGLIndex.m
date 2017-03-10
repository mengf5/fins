function [pts,lPts] = getBCGLIndex(fS,axis,side,pos)

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

bcx = bcxStart:bcxEnd;
bcy = bcyStart:bcyEnd;

bcx = (axis==0)*(bcx - (-1)^side * pos );
bcy = (axis==1)*(bcy - (-1)^side * pos );

pts = kron((bcx-1)*Nyg,ones(1,lBcx)) + kron(ones(1,lBcy),bcy);
lPts= length(pts);


end