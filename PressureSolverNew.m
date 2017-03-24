function P = PressureSolverNew(t2,fS)
%--------------------------------------------------------------------------
% sub functions contained in this code:
%
%--------------------------------------------------------------------------

%read in fS;
%--------------------------------------------------------------------------

beta  = fS.beta;
alpha = fS.alpha;

g    = fS.g;
mubx = fS.mubx;
muby = fS.muby;

x    = fS.x;
y    = fS.y;
Nxg    = fS.Nxg;
Nyg    = fS.Nyg;
hx   = fS.hx;
hy   = fS.hy;

axisX = 0;
axisY = 1;

tw = fS.tw;

UN   = fS.UN;
VN   = fS.VN;
TemN = fS.TemN;


p  = fS.p;
% fx = fS.fx;
% fy = fS.fy;
% dubcdt = fS.dubcdt;
% dvbcdt = fS.dvbcdt;
dfxdx  = fS.dfxdx;
dfydy  = fS.dfydy;

Lpl = fS.Lpl;
Lpu = fS.Lpu;
Pp  = fS.Pp;


P    = zeros(Nxg,Nyg);
RHSp = zeros(Nxg,Nyg);
%--------------------------------------------------------------------------
% fill in interior (with the boundary lines)
%--------------------------------------------------------------------------
iEnd   = ib;
jEnd   = jb;
iStart = ia;
jStart = ja;

% For periodic boundary condition,
% only the left/bottom boundary line is considered as interior

for axis =0:1
    side = 1;
    localBC = BC(axis+1,side);

    if localBC == 2
        iEnd = iEnd - (axis==0);
        jEnd = jEnd - (axis==1);
    end
end

i = iStart:iEnd;
j = jStart:jEnd;


cdx = alpha*(12*mubx);

if cdx == 0
    cdx = alpha/dtn;%1*0.1/hx^2;
end

if tw == 1
    
    forcingTerms =  dfxdx(x(i,j),y(i,j),t2) + dfydy(x(i,j),y(i,j),t2);

elseif tw == 0

    forcingTerms = 0;

end


divergence = getUx(UN,i,j,axisX,hx) + getUx(VN,i,j,axisY,hx);

RHSp(i,j) = - getUxx(UN,i,j,axisX,hx).^2 ...
	    - getUx(UN,i,j,axisY,hy).* getUx(VN,i,j,axisX,hx) ...
	    - getUxx(VN,i,j,axisY,hy).^2 ...
	    + forcingTerms ...
	    + beta*g*getUx(TemN,i,j,axisY,hy) ...
	    + cdx*divergence;

%--------------------------------------------------------------------------
% fill in ghost lines
%--------------------------------------------------------------------------

for axis = 0:1
  for side = 0:1
    
    switch localBC
                
      case 1
        % px = -ut - uux - vuy + mu (uxx+uyy) + fx
        pos = 1;
                    
        [bcx,bcy] = getBCGLlocation(axis,side,ia,ib,ja,jb,pos);
                    
        RHSp(bcx,bcy) = -getCompPx(V,fS,bcx,bcy,axis,side,pos);
          
        pos = 2;
                    
        [bcx,bcy] = getBCGLlocation(axis,side,ia,ib,ja,jb,pos);
                    
        RHSp(bcx,bcy)  = getZero(bcx,bcy);
                    
      case 3
        % p = p(t)
        for pos = 1:2;
          
          [bcx,bcy] = getBCGLlocation(axis,side,ia,ib,ja,jb,pos);
          
          RHSp(bcx,bcy) = p(x(bcx,bcy),y(bcx,bcy),t2);

        end
        
      case 4
        
      case 5
        
    end
  end
end

%--------------------------------------------------------------------------
% now fix em corner points
%--------------------------------------------------------------------------

for sideX = 0:1
  for sideY = 0:1
    
    px = (sideX==0)*ia + (sideX==1)*ib;
    py = (sideY==0)*ja + (sideY==1)*jb;
    
    px1 = px - (-1)^sideX * 1;
    py1 = py - (-1)^sideY * 1;
    
    px2 = px - (-1)^sideX * 2;
    py2 = py - (-1)^sideY * 1;
    
    px3 = px - (-1)^sideX * 1;
    py3 = py - (-1)^sideY * 2;
    
    px4 = px - (-1)^sideX * 2;
    py4 = py - (-1)^sideY * 2;


    chooseBC(1) = BC(1,1+sideX);
    chooseBC(2) = BC(2,1+sideY);

    if max( chooseBC ) < 4

      localBC=  max( chooseBC );
      
    else 

      fprintf('these could be wrong, since the developer havenot investigated this yet.\n');
      fprintf('proceed with extreme caution \n');
      pause;

    end

    switch localBC
%      case 1
% extrapolation 
%      case 2
%periodic
      case 3

	RHSp(px1,py1) = p(x(px1,py1),y(px1,py1),t2);
        
        RHSp(px2,py2) = p(x(px2,py2),y(px2,py2),t2);
        
        RHSp(px3,py3) = p(x(px3,py3),y(px3,py3),t2);
        
        RHSp(px4,py4) = p(x(px4,py4),y(px4,py4),t2);
	
	
    end	   
    
  end
end

%------------------------------------------------------------------
%from RHS vector
%------------------------------------------------------------------

RHSimp = RHSp(1,1:Nyg);
for i = 2:Nxg
  RHSimp = [RHSimp RHSp(i,1:Nyg)];
end

flagAddRow = 1;

for axis = 0:1
  for side = 0:1

    localBC = BC(axis+1,side+1);

    if localBC > 2
      flagAddRow = 0;
    end
    
  end
end

if flagAddRow == 1
 
  RHSimp = [RHSimp 0];

end

if fS.directSolve == 0 
    yt = Lpl\(Pp*RHSimp');
    P2n = Lpu\yt;

elseif fS.directSolve == 1
    Lp = fS.Lp;
    P2n = Lp\(RHSimp');
end


end




function approxUx = getCompPx(U,V,fS,bcx,bcy,axis,side,pos)

tw = fS.tw;
mu = fS.mu;
hx = fS.hx;
hy = fS.hy;

% get the index for the corresponding boundary points

iB = (axis==0)*(bcx + (-1)^side*pos) + (axis==1)*bcx;
jB = (axis==1)*(bcy + (-1)^side*pos) + (axis==0)*bcy;

% if twilightZone forcing is given, use exact dudx,
% otherwise, use approximated dudx

approxPx = 0;

if tw > 0

    u = fS.u    
    v = fS.v    
  
  if axis == 0
    dudt  =   fS.dudt;    
    dudx  = - fS.dvdy;
    dudy  =   fS.dudy;    

  elseif axis == 1
    dudt  = fS.dvdt;    
    dudx  = - fS.dudx;
    dudy  =   fS.dvdy;    

    beta = fS.beta;
    g    = fS.beta;
    TemN = fS.TemN;
    tref = fS.ref;
    
    approxPx =  beta*g*(TemN(i,j)-tref);

    

  end
  
    approxPx = approxPx - dudt(x(iB,jB),y(iB,jB),t2) ...
	       + u(x(iB,jB),y(iB,jB),t2).*dudx(x(iB,jB),y(iB,jB),t2) ...
	       - v(x(iB,jB),y(iB,jB),t2).*dudy(x(iB,jB),y(iB,jB),t2) ;

end

  if axis == 0

    f = fS.fx;
    
  elseif axis == 1

    f = fS.fy;

  end

approxPx = approxPx + getCurlCurl(U,V,iB,jB,axis,hx,hy) + f(x(iB,jB),y(iB,jB),t2);
    
end


function approxCurlCurl = getCurlCurl(U,V,iB,jB,axis,hx,hy);

  uCross = U*(axis==1) + V*(axis==0);
  uD     = U*(axis==0) + V*(axis==1);
  
  getUxy(uCross,iB,jB,hx,hy,axis) + getUxx(U,iB,jB,axis,h)
  
end


function approxUxx = getUxx(U,iB,jB,axis,h)

xShift =   (axis==0);

yShift =   (axis==1);

approxUxx = ((-1)* U(iB + xShift*2, jB + yShift*2) ...
             + 16* U(iB + xShift, jB + yShift) +   ...
             (-1)* U(iB - xShift*2, jB - yShift*2) ...
             + 16* U(iB - xShift, jB - yShift) +   ...
            (-30)* U(iB           , jB          ) )/h^2;

end

function approxUx = getUx(U,iB,jB,axis,h)

xShift =   (axis==0);

yShift =   (axis==1);

approxUx = ((-1)* U(iB + xShift*2, jB + yShift*2) ...
            + 16* U(iB + xShift, jB + yShift) + ...
            (-1)* U(iB - xShift*2, jB - yShift*2) ...
            + 16* U(iB - xShift, jB - yShift)  ...
            )/(12*h);
end


function approxUxy = getUxy(U,iB,jB,hx,hy,axis)


approxUxy = -(-U(iB+2,jB+2) + 8*U(iB+1,jB+2) - 8*U(iB-1,jB+2) + U(iB-2,jB+2)) ...
            +8*(-U(iB+2,jB+1) + 8*U(iB+1,jB+1) - 8*U(iB-1,jB+1) + U(iB-2,jB+1)) ...
            -8*(-U(iB+2,jB-1) + 8*U(iB+1,jB-1) - 8*U(iB-1,jB-1) + U(iB-2,jB-1)) ...
            + (-U(iB+2,jB-2) + 8*U(iB+1,jB-2) - 8*U(iB-1,jB-2) + U(iB-2,jB-2));

approxUxy = approxUxy/(12*hx*12*hy);

end


function zero = getZero(bcx,bcy)

zero = zeros(length(bcx),length(bcy));

end

function [bcx,bcy] = getBCGLlocation(axis,side,ia,ib,ja,jb,pos)

bcxStart = (side==0)*ia + (side==1)*ib;
bcyStart = (side==0)*ja + (side==1)*jb;

bcxEnd = (axis==1)*((side==1)*ia + (side==0)*ib) +  (axis==0)*(bcxStart);
bcyEnd = (axis==0)*((side==1)*ja + (side==0)*jb) +  (axis==1)*(bcyStart);

if bcxStart>bcxEnd
    temp = bcxEnd;
    bcxEnd = bcxStart;
    bcxStart = temp;
end

if bcyStart>bcyEnd
    temp = bcyEnd;
    bcyEnd = bcyStart;
    bcyStart = temp;
end


bcx = min(bcxStart,bcxEnd):max(bcxStart,bcxEnd);
bcy = min(bcyStart,bcyEnd):max(bcyStart,bcyEnd);
% shift from boundary to ghost line wrt its pos(ition)
bcx = bcx + (axis==0)*( - (-1)^side * pos );
bcy = bcy + (axis==1)*( - (-1)^side * pos );

end
