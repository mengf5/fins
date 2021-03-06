function P = PressureSolverNew(t2,count,fS)
%--------------------------------------------------------------------------
% sub functions contained in this code:
%
%--------------------------------------------------------------------------

%read in fS;
%--------------------------------------------------------------------------

beta  = fS.beta;
alpha = fS.alpha;

BC = fS.BC;
%BC =[3,3;3,3];

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


ia = fS.ia;
ib = fS.ib;
ja = fS.ja;
jb = fS.jb;

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
    cdx = alpha/fS.dt;%1*0.1/hx^2;
end

if tw == 1
    
    forcingTerms =  dfxdx(x(i,j),y(i,j),t2) + dfydy(x(i,j),y(i,j),t2);

elseif tw == 0

    forcingTerms = 0;

end


divergence = getUx(UN,i,j,axisX,hx) + getUx(VN,i,j,axisY,hy);

RHSp(i,j) = - getUx(UN,i,j,axisX,hx).^2 ...
	    - 2*getUx(UN,i,j,axisY,hy).* getUx(VN,i,j,axisX,hx) ...
	    - getUx(VN,i,j,axisY,hy).^2 ...
	    + forcingTerms ...
	    + beta*g*getUx(TemN,i,j,axisY,hy) ...
	    + cdx*divergence;

%RHSp(i,j) = fS.dpdx2(x(i,j),y(i,j),t2)+fS.dpdy2(x(i,j),y(i,j),t2);
%--------------------------------------------------------------------------
% fill in ghost lines
%--------------------------------------------------------------------------

for side = 0:1
  for axis = 0:1
      
      localBC = BC(axis+1,side+1);
      switch localBC
          
          case 1
              % px = -ut - uux - vuy + mu (uxx+uyy) + fx
              pos = 1;
              
              [bcx,bcy] = getBCGLlocation(axis,side,ia,ib,ja,jb,pos);
              
              RHSp(bcx,bcy) = getCompPx(UN,VN,fS,bcx,bcy,axis,side,pos,count);
        
        
%         iB = (axis==0)*(bcx + (-1)^side*pos) + (axis==1)*bcx;
%         jB = (axis==1)*(bcy + (-1)^side*pos) + (axis==0)*bcy;
% 
%         if axis == 0
%             
%             RHSp(bcx,bcy) = fS.dpdx(x(iB,jB),y(iB,jB),t2);
%             
%         elseif axis==1
%             
%             RHSp(bcx,bcy) = fS.dpdy(x(iB,jB),y(iB,jB),t2);
%             
%         end
          
        pos = 2;
                    
        [bcx,bcy] = getBCGLlocation(axis,side,ia,ib,ja,jb,pos);
                    
        RHSp(bcx,bcy)  = getZero(bcx,bcy);
%RHSp(bcx,bcy) = p(x(bcx,bcy),y(bcx,bcy),t2);
                    
      case 3
        % p = p(t)
        for pos = 0:2;
          
          [bcx,bcy] = getBCGLlocation(axis,side,ia,ib,ja,jb,pos);
          
          RHSp(bcx,bcy) = p(x(bcx,bcy),y(bcx,bcy),t2);

        end
        
      case 4
        
      case 5
        
    end
  end
end

    %----------------------------------------------------------------
    % on the physical corner and its extension (ghost points)
    %----------------------------------------------------------------
    
    for sideX = 0:1
        for sideY = 0:1
            
            pxC = (sideX==0)*ia + (sideX==1)*ib;
            pyC = (sideY==0)*ja + (sideY==1)*jb;
            
            chooseBC(1) = BC(1,1+sideX);
            chooseBC(2) = BC(2,1+sideY);
            
            if max( chooseBC ) < 4
                
                [localBC,whichBC] =  max( chooseBC );
                axisBC = (whichBC -1);
                
            else
                fprintf('these could be wrong, since the developer havenot investigated this yet.\n');
                fprintf('proceed with extreme caution \n');
                pause;
                
            end
            
            switch localBC
                
                case 1
                    %do nothing
                case 2
                    
                    side = (axisBC==0)*sideX + (axisBC==1)*sideY;
                    sideMatch = abs(side-1);
                    
                    if side == 0
                        %do nothing
                        
                    elseif side == 1
                        
                        % * = (Nx,Ny)
                        %
                        %       #
                        %       #
                        % xxxxxx*oo
                        %       x
                        %       x
                        %       x
                        
                        for pos = 0:2 % fixing *oo
                            
                            px = pxC + (axisBC == 0) * pos;
                            py = pyC + (axisBC == 1) * pos;
                            
                            RHSp(px,py) = 0;
                            
                        end
                        
                        
                        for pos = 1:2 % fixing ##
                            
                            px = pxC + (axisBC == 1) * (-1)^(sideX+1)* pos;
                            py = pyC + (axisBC == 0) * (-1)^(sideY+1)* pos;
                            
                            RHSp(px,py) = 0;
                            
                        end
                        
                        
                    end
                    
                    
                case 3
                    %do nothing
                    
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
        
        case 1
%             
            RHSp(px1,py1) = p(x(px1,py1),y(px1,py1),t2);
            
            RHSp(px2,py2) = p(x(px2,py2),y(px2,py2),t2);
            
            RHSp(px3,py3) = p(x(px3,py3),y(px3,py3),t2);
            
            RHSp(px4,py4) = p(x(px4,py4),y(px4,py4),t2);
            
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
 
    totP = sum(sum(p(x,y,t2)));
  RHSimp = [RHSimp totP];

end

if fS.directSolve == 0 
    yt = Lpl\(Pp*RHSimp');
    P2n = Lpu\yt;

elseif fS.directSolve == 1
    Lp = fS.Lp;
    P2n = Lp\(RHSimp');
end

for i = 1:Nxg
    P(i,:) = P2n((i-1)*Nyg+1:i*Nyg);
end

end




function approxPx = getCompPx(U,V,fS,bcx,bcy,axis,side,pos,count)

tw = fS.tw;
mu = fS.mu;
hx = fS.hx;
hy = fS.hy;
x  = fS.x;
y  = fS.y;
t2 = fS.t2;

% get the index for the corresponding boundary points

iB = (axis==0)*(bcx + (-1)^side*pos) + (axis==1)*bcx;
jB = (axis==1)*(bcy + (-1)^side*pos) + (axis==0)*bcy;

% if twilightZone forcing is given, use exact dudx,
% otherwise, use approximated dudx

sizePx = (axis==0)*length(jB) + (axis==1)*length(iB) ;
approxPx = zeros(1,sizePx);

if tw > 0

    u = fS.u;    
    v = fS.v;   
  
  if axis == 0
     dudt  =   fS.dudt;    
     dudx  =   fS.dvdy;
     dudy  =   fS.dudy; 
     dudx2  =   fS.dudx2; 
     dudy2  =   fS.dudy2; 

  elseif axis == 1
     dudt  =   fS.dvdt;    
     dudx  =   fS.dvdx;
     dudy  =   fS.dudx; 
     dudx2  =   fS.dvdx2; 
     dudy2  =   fS.dvdy2; 

    beta = fS.beta;
    g    = fS.g;
    TemN = fS.TemN;
    tref = fS.tref;
    
    approxPx =  beta*g*(TemN(iB,jB)-tref);

    

  end
  
    approxPx = approxPx - dudt(x(iB,jB),y(iB,jB),t2) ...
	       - (-1)^(axis + 1)*u(x(iB,jB),y(iB,jB),t2).*dudx(x(iB,jB),y(iB,jB),t2) ...
	       - (-1)^(axis    )*v(x(iB,jB),y(iB,jB),t2).*dudy(x(iB,jB),y(iB,jB),t2) ;

end

  if axis == 0

    f = fS.fx;
    
  elseif axis == 1

    f = fS.fy;

  end

% approxPx = approxPx + mu*(dudx2(x(iB,jB),y(iB,jB),t2)+dudy2(x(iB,jB),y(iB,jB),t2)) + f(x(iB,jB),y(iB,jB),t2);

% h     = (axis==0)*hx + (axis==1)*hy; 
% hP    = (axis==0)*hy + (axis==1)*hx;
% axisP = abs(axis-1);
% 
% uD = (axis==0)*U + (axis==1)*V; 
% 
% approxPx = approxPx + mu*getCurlCurl(U,V,iB,jB,axis,hx,hy) + f(x(iB,jB),y(iB,jB),t2);
% 
% if length(iB) == 1
% 
%     fixIndexI = iB;
%     fixIndexJ = [jB(1),jB(end)];
%     
% else
%     
%     fixIndexI = [iB(1),iB(end)];
%     fixIndexJ = jB;
% 
% end
% 
% indexPx = [1,length(approxPx)];

% for n = 1:2
% 
%     if length(iB) == 1
%         fixI = fixIndexI;
%         fixJ = fixIndexJ(n);
%         
%     else
%         fixI = fixIndexI(n);
%         fixJ = fixIndexJ;
%     end
%     
%     if tw>0
%     
%         approxDu = dudx2(x(fixI,fixJ),y(fixI,fixJ),t2)+dudy2(x(fixI,fixJ),y(fixI,fixJ),t2);
%     else
%         
%         approxDu = getUxx(uD,fixI,fixJ,axis,h) + getUxx(uD,fixI,fixJ,axisP,hP);
%    
%     end
% 
%     
%     approxPx(indexPx(n))   =  approxPx(indexPx(n)) ...
%         + mu*(approxDu) + f(x(fixI,fixJ),y(fixI,fixJ),t2);
%     
% end

% plot(approxPx);
% hold on
%approxPx = 0;

twilightZone = fS.twilightZone;

if length(iB) == 1
    
    if twilightZone<=6
        approxPx(1)   =  approxPx(1) + mu*(dudx2(x(iB,jB(1)),y(iB,jB(1)),t2)+dudy2(x(iB,jB(1)),y(iB,jB(1)),t2)) + f(x(iB,jB(1)),y(iB,jB(1)),t2);
    else
        approxPx(1)   =  0;
    end
    
    for i = 2:length(jB)-1
        if twilightZone<=6
            F = f(x(iB,jB(i)),y(iB,jB(i)),t2);
        else
            F = 0;
        end
        approxPx(i) = approxPx(i) + mu*getCurlCurl(U,V,iB,jB(i),axis,hx,hy) + F;
    end
    
    if twilightZone<=6
        
        approxPx(end) =  approxPx(end) + mu*(dudx2(x(iB,jB(end)),y(iB,jB(end)),t2)+dudy2(x(iB,jB(end)),y(iB,jB(end)),t2)) + f(x(iB,jB(end)),y(iB,jB(end)),t2);
    else
        approxPx(end) =  0;
    end
else
    if twilightZone<=6
        
        approxPx(1)   =  approxPx(1) + mu*(dudx2(x(iB(1),jB),y(iB(1),jB),t2)+dudy2(x(iB(1),jB),y(iB(1),jB),t2)) + f(x(iB(1),jB),y(iB(1),jB),t2);
        
    else
        approxPx(1)   =  0;
    end
    
    for i = 2:length(iB)-1
        if twilightZone<=6
            F = f(x(iB(i),jB),y(iB(i),jB),t2);
        else
            F = 0;
        end
        approxPx(i) = approxPx(i) + mu*getCurlCurl(U,V,iB(i),jB,axis,hx,hy) + F;
    end

    if twilightZone<=6
        approxPx(end) =  approxPx(end) + mu*(dudx2(x(iB(end),jB),y(iB(end),jB),t2)+dudy2(x(iB(end),jB),y(iB(end),jB),t2)) + f(x(iB(end  ),jB),y(iB(end  ),jB),t2);
    else
        approxPx(end) =  0;
    end
    
end



if axis == 0
    dpdx  =   fS.dpdx;
    
elseif axis == 1
    dpdx  =   fS.dpdy;
    
end

%approxPxE = dpdx(x(iB,jB),y(iB,jB),t2);
% plot(approxPxE - approxPx)

end


function approxCC = getCurlCurl(U,V,iB,jB,axis,hx,hy)

  uCross = U*(axis==1) + V*(axis==0);
  uD     = U*(axis==0) + V*(axis==1);
  
  axisD = abs(axis-1);
  
  approxCC = (-1)*(getUxy(uCross,iB,jB,hx,hy) - getUxx(uD,iB,jB,axisD,hy));
  
end


function approxUxx = getUxx(U,iB,jB,axis,h)

xShift =   (axis==0);

yShift =   (axis==1);

approxUxx = ((-1)* U(iB + xShift*2, jB + yShift*2) ...
             + 16* U(iB + xShift, jB + yShift) +   ...
             (-1)* U(iB - xShift*2, jB - yShift*2) ...
             + 16* U(iB - xShift, jB - yShift) +   ...
            (-30)* U(iB           , jB          ) )/(12*h^2);

end

function approxUx = getUx(U,iB,jB,axis,h)

xShift =   (axis==0);

yShift =   (axis==1);

approxUx = ((-1)* U(iB + xShift*2, jB + yShift*2) ...
            + 8 * U(iB + xShift, jB + yShift) + ...
            (+1)* U(iB - xShift*2, jB - yShift*2) ...
             -8 * U(iB - xShift, jB - yShift)  ...
            )/(12*h);
end


function approxUxy = getUxy(U,iB,jB,hx,hy)


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
