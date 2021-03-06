function [count,U,V,rhsuC,rhsvC,maxgrad] = VelocitySolverNew(t1,t2,count,fS)
%--------------------------------------------------------------------------
% sub functions contained in this code:
%1.           [UF,VF] = getFaceValues(fS,U,V,hx,hy)
%2.                 F = WENOFlux(f,fd,dir,h,i,j,Ud)
%3.  [wp1,wp2,wm1,wm2]= weight(f,dir,i,j)
%4.                 W = solveIm(Nxg,Nyg,U,dir,directSolve)


% 1.get the face value at half points like i+1/2,j
% 2.get dirivatives with BWENO scheme
% 3.get the weights used in WENOFlux
% 4.convert rhs matrix to b then solve Ax\b
%
%--------------------------------------------------------------------------

%read in fS;
%--------------------------------------------------------------------------

% if 1, use 4/6th order extrapolation on the first ghost line
% for tangetial velocity ( tangentExtOrder determines the order)
% 042317 fm
tangentExt = fS.tangentExt;

dt  = fS.dt;
% dte = fS.dte;
% dtn = fS.dtn;

tOrder    = fS.tOrder;
tMethod   = fS.tMethod;
tExplicit = fS. tExplicit;
cBI    = fS.cBI;
directSolve = fS.directSolve;
%twilightZone = fS.twilightZone;
tw = fS.tw;

beta = fS.beta;
tref = fS.tref;
g    = fS.g;
%mu   = fS.mu;
mubx = fS.mubx;
muby = fS.muby;

x    = fS.x;
y    = fS.y;
Nxg    = fS.Nxg;
Nyg    = fS.Nyg;
Nx    = fS.Nx;
%Ny    = fS.Ny;
hx   = fS.hx;
hy   = fS.hy;

BC   = fS.BC;
ad41 = fS.ad41;
ad42 = fS.ad42;


cA = fS.cA;

timeint = fS.cBE;
if count>1
    timeint = zeros(1,length(fS.cBE));
    timeint(1) = fS.cBI;
end

% bcp  = fS.bcp;
% ad21 = fS.ad21;
% ad22 = fS.ad22;
WENO = fS.WENO;
uw   = fS.uw;
cons = fS.cons;

UN = fS.UN;
UC = fS.UC;
UP1 = fS.UP1;
UP2 = fS.UP2;
UP3 = fS.UP3;

VN = fS.VN;
VC = fS.VC;
VP1 = fS.VP1;
VP2 = fS.VP2;
VP3 = fS.VP3;

PN = fS.PN;
% PC = fS.PC;
% PP1 = fS.PP1;
% PP2 = fS.PP2;
% PP3 = fS.PP3;
% PP4 = fS.Pextr;

TemN = fS.TemN;

rhsuP1=fS.rhsuP1;
rhsuP2=fS.rhsuP2;
rhsuP3=fS.rhsuP3;

rhsvP1=fS.rhsvP1;
rhsvP2=fS.rhsvP2;
rhsvP3=fS.rhsvP3;

u   =fS.u;
v   =fS.v;
% fx  =fS.fx;
% fy  =fS.fy;
% dvdt=fS.dvdt;
% dudt=fS.dudt;
% tem =fS.tem;
% dvdx=fS.dvdx;


fxE  =fS.fxE;
fyE  =fS.fyE;
fxI  =fS.fxI;
fyI  =fS.fyI;

ia = fS.ia;
ib = fS.ib;
ja = fS.ja;
jb = fS.jb;
%--------------------------------------------------------------------------

% initialization
%--------------------------------------------------------------------------
U    = zeros(Nxg,Nyg);
V    = zeros(Nxg,Nyg);
rhsU = zeros(Nxg,Nyg);
rhsV = zeros(Nxg,Nyg);

d4u1 = zeros(Nxg,Nyg);
d4v1 = zeros(Nxg,Nyg);
maxgrad = 0;
%--------------------------------------------------------------------------


% fill in interior
%--------------------------------------------------------------------------
iEnd   = ib-1;
jEnd   = jb-1;
iStart = ia+1;
jStart = ja+1;

% For periodic boundary condition,
%one boundary line is considered as interior

for axis =0:1
    side = 1;
    localBC = BC(axis+1,side);
    % on the boundary
    if localBC == 2
        
        iStart = iStart - (axis==0);
        jStart = jStart - (axis==1);
              
    end
end

% add slip wall 042417 fm
for axis = 0:1
    for side = 0:1
        
        localBC = BC(axis+1,side+1);
        
        if localBC == 4
            
            iStart = iStart - (axis==0);
            jStart = jStart - (axis==1);
            
            iEnd   = iEnd   + (axis==0);
            jEnd   = jEnd   + (axis==1);
        end
        
    end
end

i = iStart:iEnd;
j = jStart:jEnd;

% calculate artifitial dissipation
%--------------------------------------------------------------------------
if max(abs(ad41),abs(ad42)) > 0
    
    graduv1 = 0.25*(abs((-UN(i+2,j) + 8*UN(i+1,j) - 8*UN(i-1,j) + UN(i-2,j))/(12*hx)) + ...
        abs((-UN(i,j+2) + 8*UN(i,j+1) - 8*UN(i,j-1) + UN(i,j-2))/(12*hy))+...
        abs((-VN(i+2,j) + 8*VN(i+1,j) - 8*VN(i-1,j) + VN(i-2,j))/(12*hx))+...
        abs((-VN(i,j+2) + 8*VN(i,j+1) - 8*VN(i,j-1) + VN(i,j-2))/(12*hy)));
    maxgrad = max(max(graduv1));
    
    lapu1 = (UN(i+2,j) -4*UN(i+1,j) +6*UN(i,j) -4*UN(i-1,j)+UN(i-2,j)) + ...
        (UN(i,j+2) -4*UN(i,j+1) +6*UN(i,j) -4*UN(i,j-1)+UN(i,j-2));
    
    lapv1 = (VN(i+2,j) -4*VN(i+1,j) +6*VN(i,j) -4*VN(i-1,j)+VN(i-2,j)) + ...
        (VN(i,j+2) -4*VN(i,j+1) +6*VN(i,j) -4*VN(i,j-1)+VN(i,j-2));
    
    d4u1(i,j) = -(ad41+graduv1*ad42).*lapu1;
    d4v1(i,j) = -(ad41+graduv1*ad42).*lapv1;
end
%--------------------------------------------------------------------------

if cons == 1
    % non-conservative form cons == 1
    if WENO == 1
        % BWENO
        dirx = [1,0];
        diry = [0,1];
        
        [UF,VF] = getFaceValues(fS,UN,VN,hx,hy);
        
        ux = WENOFlux(UN,UN,dirx,hx,i,j,UF,uw);
        uy = WENOFlux(UN,VN,diry,hy,i,j,VF,uw);
        
        vx = WENOFlux(VN,UN,dirx,hx,i,j,UF,uw);
        vy = WENOFlux(VN,VN,diry,hy,i,j,VF,uw);
    else
        % centered difference
        ux(i,j) = (-UN(i+2,j) + 8*UN(i+1,j) - 8*UN(i-1,j) + UN(i-2,j))/(12*hx);
        uy(i,j) = (-UN(i,j+2) + 8*UN(i,j+1) - 8*UN(i,j-1) + UN(i,j-2))/(12*hy);
        vx(i,j) = (-VN(i+2,j) + 8*VN(i+1,j) - 8*VN(i-1,j) + VN(i-2,j))/(12*hx);
        vy(i,j) = (-VN(i,j+2) + 8*VN(i,j+1) - 8*VN(i,j-1) + VN(i,j-2))/(12*hy);
    end
    
    convu(i,j) = UN(i,j).*ux(i,j) + VN(i,j).*uy(i,j);
    convv(i,j) = UN(i,j).*vx(i,j) + VN(i,j).*vy(i,j);
    
else
    % divR
    div21u1 = ((1/2*(UN(i+1,j)+UN(i,j))).^2 - (1/2*(UN(i,j)+UN(i-1,j))).^2)/(hx) + ...
        ((1/2*(VN(i,j+1)+VN(i,j))).*(1/2*(UN(i,j+1)+UN(i,j))) ...
        - (1/2*(VN(i,j-1)+VN(i,j))).*(1/2*(UN(i,j-1)+UN(i,j))))/(hy);
    div22u1 = ((1/2*(UN(i+2,j)+UN(i,j))).^2 - (1/2*(UN(i,j)+UN(i-2,j))).^2)/(2*hx) + ...
        ((1/2*(VN(i,j+2)+VN(i,j))).*(1/2*(UN(i,j+2)+UN(i,j))) ...
        - (1/2*(VN(i,j-2)+VN(i,j))).*(1/2*(UN(i,j-2)+UN(i,j))))/(2*hy);
    
    convu(i,j) = 4/3*div21u1 - 1/3*div22u1;
    
    div21v1 = ((1/2*(VN(i,j+1)+VN(i,j))).^2 - (1/2*(VN(i,j)+VN(i,j-1))).^2)/(hy) + ...
        ((1/2*(VN(i+1,j)+VN(i,j))).*(1/2*(UN(i+1,j)+UN(i,j))) ...
        - (1/2*(VN(i-1,j)+VN(i,j))).*(1/2*(UN(i-1,j)+UN(i,j))))/(hx);
    
    div22v1 = ((1/2*(VN(i,j+2)+VN(i,j))).^2 - (1/2*(VN(i,j)+VN(i,j-2))).^2)/(2*hy) + ...
        ((1/2*(VN(i+2,j)+VN(i,j))).*(1/2*(UN(i+2,j)+UN(i,j))) ...
        - (1/2*(VN(i-2,j)+VN(i,j))).*(1/2*(UN(i-2,j)+UN(i,j))))/(2*hx);
    
    convv(i,j) = 4/3*div21v1 - 1/3*div22v1;
    
    if cons==3
        % skew symmetric
        cont1(i,j) = (UN(i+1,j)-UN(i-1,j))/(2*hx) + (VN(i,j+1)-VN(i,j-1))/(2*hy);
        cont2(i,j) = (UN(i+2,j)-UN(i-2,j))/(4*hx) + (VN(i,j+2)-VN(i,j-2))/(4*hy);
        
        div21u1 = div21u1 - 0.5*UN(i,j).*cont1(i,j);
        div22u1 = div22u1 - 0.5*UN(i,j).*cont2(i,j);
        
        convu(i,j) = 4/3*div21u1 - 1/3*div22u1;
        
        div21v1 = div21v1 - 0.5*VN(i,j).*cont1(i,j);
        div22v1 = div22v1 - 0.5*VN(i,j).*cont2(i,j);
        
        convv(i,j) = 4/3*div21v1 - 1/3*div22v1;
        
    end
    
end

% form the rhs matrix of the convective part fE - convU - p_x
%--------------------------------------------------------------------------
if tw == 1
    
    fxE = fxE(x(i,j),y(i,j),t1);
    fyE = fyE(x(i,j),y(i,j),t1);
    
elseif tw==0
    
    if fS.twilightZone>=7
        fxE = fxE(x(i,j),y(i,j),t1);
    else
        fxE = 0;
    end
    
    fyE = 0;
    
end


rhsuC =  fxE ...
    + d4u1(i,j) + ...
    - convu(i,j)...
    - (-PN(i+2,j) + 8*PN(i+1,j) - 8*PN(i-1,j) + PN(i-2,j))/(12*hx);

rhsvC = beta*g*(TemN(i,j)-tref) + fyE ...
    + d4v1(i,j) + ...
    - convv(i,j) ...
    - (-PN(i,j+2) + 8*PN(i,j+1) - 8*PN(i,j-1) + PN(i,j-2))/(12*hy);
%--------------------------------------------------------------------------

% form the rhs matrix of the diffusion part fI + mn*laplacian(U)
%--------------------------------------------------------------------------
if tw==1
    
    fxC = fxI(x(i,j),y(i,j),t1);
    fyC = fyI(x(i,j),y(i,j),t1);
    
elseif tw==0
    
    fxC=0;
    fyC=0;
    
end

rhsuImC = mubx*(-UN(i+2,j) + 16*UN(i+1,j) -30*UN(i,j) + 16*UN(i-1,j)-UN(i-2,j)) ...
    + muby*(-UN(i,j+2) + 16*UN(i,j+1) -30*UN(i,j) + 16*UN(i,j-1)-UN(i,j-2))...
    + fxC;

rhsvImC =  mubx*(-VN(i+2,j) + 16*VN(i+1,j) -30*VN(i,j) + 16*VN(i-1,j)-VN(i-2,j))...
    + muby*(-VN(i,j+2) + 16*VN(i,j+1) -30*VN(i,j) + 16*VN(i,j-1)-VN(i,j-2))...
    + fyC;

if tExplicit==1
    
    rhsuC = rhsuC + rhsuImC;
    
    rhsvC = rhsvC + rhsvImC;
    
end
%--------------------------------------------------------------------------
% count = 0 is the start up stage, if count =0, rhsuC, rhsvC are all we need
% count = 1 is the predictor
% count > 1 is the corrector

if count == 0
    return
end

%--------------------------------------------------------------------------

if tExplicit==1
    % explicit time stepper
    %--------------------------------------------------------------------------
    if tOrder==2
        
        rhsU(i,j) = timeint(1)*rhsuC...
            + timeint(2)*rhsuP1...
            + UC(i,j); %
        
        rhsV(i,j) = timeint(1)*rhsvC...
            + timeint(2)*rhsvP1...
            + VC(i,j);
        
        
    elseif tOrder==4
        if tMethod==1
            rhsU(i,j) = timeint(1)*rhsuC...
                + timeint(2)*rhsuP1...
                + timeint(3)*rhsuP2...
                + timeint(4)*rhsuP3...
                + UC(i,j); %
            
            rhsV(i,j) = timeint(1)*rhsvC...
                + timeint(2)*rhsvP1...
                + timeint(3)*rhsvP2...
                + timeint(4)*rhsvP3...
                + VC(i,j);
            
        elseif tMethod==2
            
            rhsU(i,j) = timeint(1)*rhsuC...
                + timeint(2)*rhsuP1...
                + timeint(3)*rhsuP2...
                + timeint(4)*rhsuP3...
                - ( cA(1)  *UC(i,j) ...
                +cA(2)  *UP1(i,j) ...
                +cA(3)  *UP2(i,j) ...
                +cA(4)  *UP3(i,j)); %
            
            
            
            rhsV(i,j) = timeint(1)*rhsvC...
                + timeint(2)*rhsvP1...
                + timeint(3)*rhsvP2...
                + timeint(4)*rhsvP3...
                - ( cA(1)  *VC(i,j) ...
                +cA(2)  *VP1(i,j) ...
                +cA(3)  *VP2(i,j) ...
                +cA(4)  *VP3(i,j));
        end
        
    end
    
    
elseif tExplicit==0
    % implicit time stepper
    %-------------------------------------------------------------------------------
    if count>1
        
        if tw == 1
            fxC = fxI(x(i,j),y(i,j),t2-dt);
            fyC = fyI(x(i,j),y(i,j),t2-dt);
        elseif tw==0
            fxC=0;
            fyC=0;
        end
        
        rhsuImC = mubx*(-UC(i+2,j) + 16*UC(i+1,j) -30*UC(i,j) + 16*UC(i-1,j)-UC(i-2,j)) ...
            + muby*(-UC(i,j+2) + 16*UC(i,j+1) -30*UC(i,j) + 16*UC(i,j-1)-UC(i,j-2))...
            + fxC;
        
        rhsvImC =  mubx*(-VC(i+2,j) + 16*VC(i+1,j) -30*VC(i,j) + 16*VC(i-1,j)-VC(i-2,j))...
            + muby*(-VC(i,j+2) + 16*VC(i,j+1) -30*VC(i,j) + 16*VC(i,j-1)-VC(i,j-2))...
            + fyC;
    end
    
    if tOrder==2
        
        if tw== 1;
            fxN = fxI(x(i,j),y(i,j),t2);
            fyN = fyI(x(i,j),y(i,j),t2);
        elseif tw==0
            fxN = 0;
            fyN = 0;
        end
        
        rhsU(i,j) = timeint(1)*rhsuC...
            + timeint(2)*rhsuP1...
            + UC(i,j)...
            + cBI(2)*rhsuImC ...
            + cBI(1)*fxN;
        
        rhsV(i,j) = timeint(1)*rhsvC...
            + timeint(2)*rhsvP1...
            + VC(i,j)...
            + cBI(2)*rhsvImC...
            + cBI(1)*fyN;
        
    elseif tOrder==4
        
        
        if tMethod==1
            
            if tw == 1
                totalFxI = cBI(1)*fxI(x(i,j),y(i,j),t2) ;
                totalFyI = cBI(1)*fyI(x(i,j),y(i,j),t2) ;
                
                fxP1 = fxI(x(i,j),y(i,j),t2-2*dt);
                fyP1 = fyI(x(i,j),y(i,j),t2-2*dt);
                
                fxP2 = fxI(x(i,j),y(i,j),t2-3*dt);
                fyP2 = fyI(x(i,j),y(i,j),t2-3*dt);
                
            elseif tw == 0
                
                totalFxI = 0;
                totalFyI = 0;
                
                fxP1 = 0;
                fyP1 = 0;
                
                fxP2 = 0;
                fyP2 = 0;
                
            end
            
            rhsuImP1 = mubx*(-UP1(i+2,j) + 16*UP1(i+1,j) -30*UP1(i,j) + 16*UP1(i-1,j)-UP1(i-2,j)) ...
                + muby*(-UP1(i,j+2) + 16*UP1(i,j+1) -30*UP1(i,j) + 16*UP1(i,j-1)-UP1(i,j-2))...
                + fxP1;
            
            rhsvImP1 =  mubx*(-VP1(i+2,j) + 16*VP1(i+1,j) -30*VP1(i,j) + 16*VP1(i-1,j)-VP1(i-2,j))...
                + muby*(-VP1(i,j+2) + 16*VP1(i,j+1) -30*VP1(i,j) + 16*VP1(i,j-1)-VP1(i,j-2))...
                + fyP1;
            
            rhsuImP2 = mubx*(-UP2(i+2,j) + 16*UP2(i+1,j) -30*UP2(i,j) + 16*UP2(i-1,j)-UP2(i-2,j)) ...
                + muby*(-UP2(i,j+2) + 16*UP2(i,j+1) -30*UP2(i,j) + 16*UP2(i,j-1)-UP2(i,j-2))...
                + fxP2;
            
            rhsvImP2 =  mubx*(-VP2(i+2,j) + 16*VP2(i+1,j) -30*VP2(i,j) + 16*VP2(i-1,j)-VP2(i-2,j))...
                + muby*(-VP2(i,j+2) + 16*VP2(i,j+1) -30*VP2(i,j) + 16*VP2(i,j-1)-VP2(i,j-2))...
                + fyP2;
            
            
            rhsU(i,j) = timeint(1)*rhsuC...
                + timeint(2)*rhsuP1...
                + timeint(3)*rhsuP2...
                + timeint(4)*rhsuP3...
                + UC(i,j)...
                + cBI(2)*rhsuImC...
                + cBI(3)*rhsuImP1...
                + cBI(4)*rhsuImP2...
                + totalFxI; %
            
            rhsV(i,j) = timeint(1)*rhsvC...
                + timeint(2)*rhsvP1...
                + timeint(3)*rhsvP2...
                + timeint(4)*rhsvP3...
                + VC(i,j)...
                + cBI(2)*rhsvImC...
                + cBI(3)*rhsvImP1...
                + cBI(4)*rhsvImP2...
                + totalFyI;
            
        elseif tMethod==2
            
            if tw == 1
                
                totalFxI = cBI(1)*fxI(x(i,j),y(i,j),t2);
                totalFyI = cBI(1)*fyI(x(i,j),y(i,j),t2);
                
            elseif tw==0
                
                totalFxI = 0;
                totalFyI = 0;
                
            end
            
            rhsU(i,j) = timeint(1)*rhsuC...
                + timeint(2)*rhsuP1...
                + timeint(3)*rhsuP2...
                + timeint(4)*rhsuP3...
                - ( cA(1)  *UC(i,j) ...
                +cA(2)  *UP1(i,j) ...
                +cA(3)  *UP2(i,j) ...
                +cA(4)  *UP3(i,j))...
                +totalFxI; %
            
            rhsV(i,j) = timeint(1)*rhsvC...
                + timeint(2)*rhsvP1...
                + timeint(3)*rhsvP2...
                + timeint(4)*rhsvP3...
                - ( cA(1)  *VC(i,j) ...
                +cA(2)  *VP1(i,j) ...
                +cA(3)  *VP2(i,j) ...
                +cA(4)  *VP3(i,j))...
                + totalFyI;
        end
        
    end
    
end


% calculating boundary and ghost points


if  tExplicit == 1
    % explicit time stepping
    %-----------------------------------------------------------------------------------
    msg = 'boundary conditions for explicit time stepping method is under construction';
    error(msg)
    
    %     for axis = 0:1
    %         for side = 0:1
    %
    %             localBC = BC(axis+1,side+1);
    %             % on the boundary
    %             pos = 0;
    %
    %             [bcx,bcy] = getBCGLlocation(axis,side,ia,ib,ja,jb,pos);
    %
    %             switch localBC
    %                 case 1
    %
    %             end
    %
    %         end
    %     end
    
elseif tExplicit == 0
    % implicit time stepping
    %-----------------------------------------------------------------------------------
    for side = 0:1
        for axis = 0:1
            
            localBC = BC(axis+1,side+1);
            % on the boundary
            pos = 0;
            
            [bcx,bcy] = getBCGLlocation(axis,side,ia,ib,ja,jb,pos);
            
            switch localBC
                
                case 1
                    % u = u(t)
                    rhsU(bcx,bcy) = u(x(bcx,bcy),y(bcx,bcy),t2);
                    rhsV(bcx,bcy) = v(x(bcx,bcy),y(bcx,bcy),t2);
                    
                case 2
                    
                case 3
                    
                    % u = u(t)
                    rhsU(bcx,bcy) = u(x(bcx,bcy),y(bcx,bcy),t2);
                    rhsV(bcx,bcy) = v(x(bcx,bcy),y(bcx,bcy),t2);
                    
                case 4
                             
                    
                case 5
                    
            end
            %--------------------------------------------------------------
            % on the ghost lines
            %--------------------------------------------------------------
            switch localBC
                
                case 1
                    
                    % u_x = -v_y
                    pos = 1;
                    
                    [bcx,bcy] = getBCGLlocation(axis,side,ia,ib,ja,jb,pos);
                    
                    if axis == 0
                        
                        rhsU(bcx,bcy) =  getCompUx(V,fS,bcx,bcy,axis,side,pos);
                        axisV=1;
                        
                        if tangentExt == 1
                            
                            rhsV(bcx,bcy)  = getZero(bcx,bcy);
                            
                        elseif tangentExt ==0
                            
                            [bcxB,bcyB] = getBCGLlocation(axis,side,ia,ib,ja,jb,0);
                            mu4 = getArtifitialDissipation( fS,axis,bcxB,bcyB,hx,hy);
                            %mu4 = 1e-3;
                            
                            rhsV(bcx,bcy) =  getCompUxx(V,fS,bcx,bcy,axisV,side,pos,count,mu4);
                            
                        end
                        
                    elseif axis == 1
                        
                        rhsV(bcx,bcy) =  getCompUx(U,fS,bcx,bcy,axis,side,pos);
                        axisU = 0;
                        
                        if tangentExt == 1
                            
                            rhsU(bcx,bcy) =  getZero(bcx,bcy);
                            
                        elseif tangentExt == 0
                            
                            [bcxB,bcyB] = getBCGLlocation(axis,side,ia,ib,ja,jb,0);
                            mu4 = getArtifitialDissipation( fS,axis,bcxB,bcyB,hx,hy);
                            %mu4 = 1e-3;
                            
                            rhsU(bcx,bcy) =  getCompUxx(U,fS,bcx,bcy,axisU,side,pos,count,mu4);
                            
                        end
                    end
                    
                    pos = 2;
                    
                    [bcx,bcy] = getBCGLlocation(axis,side,ia,ib,ja,jb,pos);
                    
                    iB = (axis==0)*(bcx + (-1)^side*pos) + (axis==1)*bcx;
                    jB = (axis==1)*(bcy + (-1)^side*pos) + (axis==0)*bcy;
                    
                    
                    if axis == 0
                        rhsU(bcx,bcy)  = getCompUxxx(fS,count,iB,jB);
                        rhsV(bcx,bcy)  = getZero(bcx,bcy);
                        
                    elseif axis == 1
                        rhsU(bcx,bcy)  = getZero(bcx,bcy);
                        rhsV(bcx,bcy)  = -getUxy(U,iB,jB,hx,hy);
                        
                    end
                    
                case 2
                    
                    
                case 3
                    % u = u(t)
                    for pos = 1:2;
                        
                        [bcx,bcy] = getBCGLlocation(axis,side,ia,ib,ja,jb,pos);
                        
                        rhsU(bcx,bcy) = u(x(bcx,bcy),y(bcx,bcy),t2);
                        rhsV(bcx,bcy) = v(x(bcx,bcy),y(bcx,bcy),t2);
                        
                    end
                    
                case 4
                    
                    pos = 1;
                    
                    [bcx,bcy] = getBCGLlocation(axis,side,ia,ib,ja,jb,pos);
                    
                    if axis == 0
                        
                        rhsU(bcx,bcy) =  getCompUx(V,fS,bcx,bcy,axis,side,pos);

                        rhsV(bcx,bcy) =  getCompUxZero(fS,bcx,bcy,axis,side,pos);
                    
                    elseif axis == 1
                        
                        rhsU(bcx,bcy) =  getCompUxZero(fS,bcx,bcy,axis,side,pos);

                    end
                    
                    
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
                            
                            rhsU(px,py) = 0;
                            rhsV(px,py) = 0;
                            
                        end
                        
                        
                        for pos = 1:2 % fixing ##
                            
                            px = pxC + (axisBC == 1) * (-1)^(sideX+1)* pos;
                            py = pyC + (axisBC == 0) * (-1)^(sideY+1)* pos;
                            
                            rhsU(px,py) = 0;
                            rhsV(px,py) = 0;
                            
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
                
                localBC =  max( chooseBC );
                
            else
                
                fprintf('these could be wrong, since the developer havenot investigated this yet.\n');
                fprintf('proceed with extreme caution \n');
                pause;
                
            end
            
            switch localBC
                
                case 1
                    
                    rx1 =  (-1)^(sideX);
                    ry1 =  (-1)^(sideY);
                    
                    rx2 =  (-1)^(sideX)*(2);
                    ry2 =  (-1)^(sideY);
                    
                    rx3 =  (-1)^(sideX);
                    ry3 =  (-1)^(sideY)*(2);
                    
                    coeff  =[3/2 , 3];
                    coeff4 =[24  , 24];
                    
                    [D1u,D2u,D1v,D2v]  = getDs(fS,px,py,rx1,ry1,hx,hy);
                    
                    rhsU(px1,py1) = coeff(1)*D1u +  coeff(2)*D2u;
                    rhsV(px1,py1) = coeff(1)*D1v +  coeff(2)*D2v;
                    
                    rhsU(px4,py4) = coeff4(1)*D1u +  coeff4(2)*D2u;
                    rhsV(px4,py4) = coeff4(1)*D1v +  coeff4(2)*D2v;
                    
                    [D1u,D2u,D1v,D2v]  = getDs(fS,px,py,rx2,ry2,hx,hy);
                    
                    rhsU(px2,py2) = coeff(1)*D1u +  coeff(2)*D2u;
                    rhsV(px2,py2) = coeff(1)*D1v +  coeff(2)*D2v;
                    
                    [D1u,D2u,D1v,D2v]  = getDs(fS,px,py,rx3,ry3,hx,hy);
                    
                    rhsU(px3,py3) = coeff(1)*D1u +  coeff(2)*D2u;
                    rhsV(px3,py3) = coeff(1)*D1v +  coeff(2)*D2v;
                    
                    %                     rhsU(px1,py1) = u(x(px1,py1),y(px1,py1),t2);
                    %                     rhsV(px1,py1) = v(x(px1,py1),y(px1,py1),t2);
                    %
                    %                     rhsU(px2,py2) = u(x(px2,py2),y(px2,py2),t2);
                    %                     rhsV(px2,py2) = v(x(px2,py2),y(px2,py2),t2);
                    %
                    %                     rhsU(px3,py3) = u(x(px3,py3),y(px3,py3),t2);
                    %                     rhsV(px3,py3) = v(x(px3,py3),y(px3,py3),t2);
                    %
                    %                     rhsU(px4,py4) = u(x(px4,py4),y(px4,py4),t2);
                    %                     rhsV(px4,py4) = v(x(px4,py4),y(px4,py4),t2);
                    
                    
                case 3
                    
                    rhsU(px1,py1) = u(x(px1,py1),y(px1,py1),t2);
                    rhsV(px1,py1) = v(x(px1,py1),y(px1,py1),t2);
                    
                    rhsU(px2,py2) = u(x(px2,py2),y(px2,py2),t2);
                    rhsV(px2,py2) = v(x(px2,py2),y(px2,py2),t2);
                    
                    rhsU(px3,py3) = u(x(px3,py3),y(px3,py3),t2);
                    rhsV(px3,py3) = v(x(px3,py3),y(px3,py3),t2);
                    
                    rhsU(px4,py4) = u(x(px4,py4),y(px4,py4),t2);
                    rhsV(px4,py4) = v(x(px4,py4),y(px4,py4),t2);
                    
                    
            end
            
            
            
            %             px1 = [px1,px1E];
            %             px2 = [px2,px2E];
            %             px3 = [px3,px3E];
            %             px4 = [px4,px4E];
            %
            %             py1 = [py1,py1E];
            %             py2 = [py2,py2E];
            %             py3 = [py3,py3E];
            %             py4 = [py4,py4E];
            %
            %             for i = 1:length(px1)
            %
            %             end
            
        end
    end
    
    % solve the problem
    %-----------------------------------------------------------------------------
    axis = 0;
    U = solveIm(fS,Nxg,Nyg,rhsU,axis,directSolve);
    
    %fix the second ghost line for V for BC 1
    %----------------------------------------------------------------------------
    
    axis = 1;
    for side = 0:1
        switch localBC
            case 1
                
                pos = 2;
                
                [bcx,bcy] = getBCGLlocation(axis,side,ia,ib,ja,jb,pos);
                iB = (axis==0)*(bcx + (-1)^side*pos) + (axis==1)*bcx;
                jB = (axis==1)*(bcy + (-1)^side*pos) + (axis==0)*bcy;
                
                rhsV(bcx,bcy)  = -getUxy(U,iB,jB,hx,hy);
                
        end
    end
    
    
    % solve the problem
    %-----------------------------------------------------------------------------
    axis = 1;
    V = solveIm(fS,Nxg,Nyg,rhsV,axis,directSolve);
    
    
    
end

count = count + 1;

end
function [UF,VF] = getFaceValues(fS,U,V,hx,hy)

Nxg =fS.Nxg;
Nyg =fS.Nyg;

for m = 3:Nxg-2
    for n = 3:Nyg-2
        
        UF(m,n) = U(m,n) ...
            + (hx/2) * (-U(m+2,n) + 8*U(m+1,n) - 8*U(m-1,n) + U(m-2,n))/(12*hx);
        VF(m,n) = V(m,n) ...
            + (hy/2) * (-V(m,n+2) + 8*V(m,n+1) - 8*V(m,n-1) + V(m,n-2))/(12*hy);
        
    end
end

m = 3;%taylor series is centered on this point//
%    m=2 point is calculated
n = 3:Nyg-2;
UF(m-1,n) = U(m,n) ...
    - (hx/2) * (-U(m+2,n) + 8*U(m+1,n) - 8*U(m-1,n) + U(m-2,n))/(12*hx);

m = 3:Nxg-2;
n = 3;
VF(m,n-1) = V(m,n) ...
    - (hy/2) * (-V(m,n+2) + 8*V(m,n+1) - 8*V(m,n-1) + V(m,n-2))/(12*hy);


end

function F = WENOFlux(f,fd,dir,h,i,j,Ud,uw)
lx = i(1):i(end);
ly = j(1):j(end);
[wplw(lx,ly),wprw(lx,ly),wmlw(lx,ly),wmrw(lx,ly)]...
    = weight(f,dir,i,j);



fp(i,j) = fd(i,j) > 0;
fm(i,j) = fd(i,j) <= 0;

fPp = fp;
fPm = fm;

fMp = fp;
fMm = fm;

for m = i(1):i(end)
    for n = j(1):j(end);
        
        if Ud(m,n) -eps > 0  %this is the +1/2 face
            fPp(m,n) = 1;
            fPm(m,n) = 0;
        end
        
        if Ud(m-dir(1),n-dir(2))-eps > 0  %this is the -1/2 face
            fMp(m,n) = 1;
            fMm(m,n) = 0;
        end
        
    end
end


if uw == 0
    wpl = max(wplw,wprw).*fPp + min(wplw,wprw).*fPm;
    wpr = max(wplw,wprw).*fPm + min(wplw,wprw).*fPp;
    
    wml = max(wmlw,wmrw).*fMp + min(wmlw,wmrw).*fMm;
    wmr = max(wmlw,wmrw).*fMm + min(wmlw,wmrw).*fMp;
else
    wpl = fp;
    wpr = fm;
    
    wml = fp;
    wmr = fm;
end

Fpl(i,j) = 1/6*( -f(i-1*dir(1),j-1*dir(2)) + 5*f(i,j)                   +2*f(i+1*dir(1),j+1*dir(2)));
Fpr(i,j) = 1/6*(2*f(i,j)                   + 5*f(i+1*dir(1),j+1*dir(2)) -  f(i+2*dir(1),j+2*dir(2)));
Fml(i,j) = 1/6*( -f(i-2*dir(1),j-2*dir(2)) + 5*f(i-1*dir(1),j-1*dir(2)) +2*f(i,j));
Fmr(i,j) = 1/6*(2*f(i-1*dir(1),j-1*dir(2)) + 5*f(i,j)                   -  f(i+1*dir(1),j+1*dir(2)));

Fp(i,j) = wpl(i,j).*Fpl(i,j) + wpr(i,j).*Fpr(i,j);
Fm(i,j) = wml(i,j).*Fml(i,j) + wmr(i,j).*Fmr(i,j);

F(i,j) = (Fp(i,j) - Fm(i,j))/h;
end

function [wp1,wp2,wm1,wm2]= weight(f,dir,i,j)

Apl   = f(i+1*dir(1),j+1*dir(2))-2*f(i,j)                  +f(i-1*dir(1),j-1*dir(2));
Bpl   = f(i+1*dir(1),j+1*dir(2))-  f(i-1*dir(1),j-1*dir(2));
Apr   = f(i+2*dir(1),j+2*dir(2))-2*f(i+1*dir(1),j+1*dir(2))+f(i,j);
Bpr   = f(i+2*dir(1),j+2*dir(2))-  f(i,j);

Aml   = f(i,j)                  -2*f(i-1*dir(1),j-1*dir(2))+f(i-2*dir(1),j-2*dir(2));
Bml   = f(i,j)                  -  f(i-2*dir(1),j-2*dir(2));
Amr   = f(i+1*dir(1),j+1*dir(2))-2*f(i,j)                  +f(i-1*dir(1),j-1*dir(2));
Bmr   = f(i+1*dir(1),j+1*dir(2))-  f(i-1*dir(1),j-1*dir(2));


betapl= 4/3*Apl.^2 + 1/2*Apl.*Bpl + 1/4*Bpl.^2;
betapr= 4/3*Apr.^2 - 1/2*Apr.*Bpr + 1/4*Bpr.^2;

betaml= 4/3*Aml.^2 + 1/2*Aml.*Bml + 1/4*Bml.^2;
betamr= 4/3*Amr.^2 - 1/2*Amr.*Bmr + 1/4*Bmr.^2;

ep    = 1e-15;
ap1    = (1/2)./(ep + betapl).^(4);
ap2    = (1/2)./(ep + betapr).^(4);
am1    = (1/2)./(ep + betaml).^(4);
am2    = (1/2)./(ep + betamr).^(4);

wp1 = ap1./(ap1+ap2);
wp2 = ap2./(ap1+ap2);

wm1 = am1./(am1+am2);
wm2 = am2./(am1+am2);

bp1 = wp1.*(3/4+wp1.*(wp1-1/2));
bp2 = wp2.*(3/4+wp2.*(wp2-1/2));
bm1 = wm1.*(3/4+wm1.*(wm1-1/2));
bm2 = wm2.*(3/4+wm2.*(wm2-1/2));

wp1 = bp1./(bp1+bp2);
wp2 = bp2./(bp1+bp2);

wm1 = bm1./(bm1+bm2);
wm2 = bm2./(bm1+bm2);

end

function W = solveIm(fS,Nxg,Nyg,U,axis,directSolve)
%% reshape to the rhs vector
rhsu =zeros(Nxg*Nyg,1);

for l = 1:Nxg
    
    index = ((l-1)*Nyg+1:l*Nyg);
    J     = 1:Nyg;
    rhsu(index) = U(l,J);
    
end

%% get the LHS matrix
if directSolve == 0
    
    if axis == 0
        Ll = fS.Lul;
        Lu = fS.Luu;
        p  = fS.Up;
    elseif axis ==1
        Ll = fS.Lvl;
        Lu = fS.Lvu;
        p  = fS.Vp;
    end
    yt = Ll\(p*rhsu);
    Ut  = Lu\yt;
    
elseif directSolve == 1
    
    if axis == 0
        
        Lu = fS.Lu;
        Ut  = Lu\rhsu;
        
    elseif axis ==1
        
        Lu = fS.Lv;
        Ut  = Lu\rhsu;
    end
    
end

%% change the vector U,V back to matrix
for l = 1:Nxg
    index = ((l-1)*Nyg+1:l*Nyg);
    W(l,:) = Ut(index);
end
end


function Uxxx = getCompUxxx(fS,count,i,j)
tw = fS.tw;
t2 = fS.t2;
x = fS.x;
y = fS.y;

%dvdx2y = fS.dvdx2y;
%Uxxx  = -dvdx2y(x(i,j),y(i,j),t2);
%return

% dudx3 = -dvdx2y
%       =  dvdy3 - 1/mu*(dvdty + dudy*dvdx + u*dvdxy
%                       +dvdy*dvdy + v*dvdy2 + dpdy2);

mu    = fS.mu;
u = fS.u;
v = fS.v;
dudy  = fS.dudy;
dvdy3 = fS.dvdy3;
dvdy2 = fS.dvdy2;
dvdx  = fS.dvdx;
dvdy  = fS.dvdy;
dvdxy = fS.dvdxy;
dvdty = fS.dvdty;
dfydy = fS.dfydy;
order=2;
hy = fS.hy;

dpdy2Approx = getDpdx(fS,i,j,count,1,hy,order);
%dpdy2Approx = fS.dpdy2(x(i,j),y(i,j),t2);

%dudx3 = -dvdx2y(x(i,j),y(i,j),t2);
if tw == 1
    
    Uxxx = dvdy3(x(i,j),y(i,j),t2) - 1/mu*(dvdty(x(i,j),y(i,j),t2) ...
        + dudy(x(i,j),y(i,j),t2).*dvdx(x(i,j),y(i,j),t2) ...
        + u(x(i,j),y(i,j),t2).*dvdxy(x(i,j),y(i,j),t2) ...
        + dvdy(x(i,j),y(i,j),t2).*dvdy(x(i,j),y(i,j),t2) ...
        + v(x(i,j),y(i,j),t2).*dvdy2(x(i,j),y(i,j),t2) ...
        + dpdy2Approx - dfydy(x(i,j),y(i,j),t2));
    
elseif tw == 0
    
    Uxxx =  - 1/mu*( ...
        + dpdy2Approx - dfydy(x(i,j),y(i,j),t2));
    
end

end


function Uxx = getCompUxx(U,fS,bcx,bcy,axis,side,pos,count,mu4)

tw = fS.tw;
mu = fS.mu;
t2 = fS.t2;

x = fS.x;
y = fS.y;

mu = mu + mu4;

% keep the exact uxx here for debug
% dudx2 = (axis==0)*fS.dudx2 + (axis=1)*fS.dvdy2;
% U(bcx,bcy) = dvdx2(x(iB,jB),y(iB,jB),t2);

% get the index for the corresponding boundary points

iB = (axis==1)*(bcx + (-1)^side*pos) + (axis==0)*bcx;
jB = (axis==0)*(bcy + (-1)^side*pos) + (axis==1)*bcy;

hx = fS.hx;
hy = fS.hy;
h  = (axis==0)*hx  + (axis==1)*hy;

order = 1;
dpdxApprox = getDpdx(fS,iB,jB,count,axis,h,order);



% get the forcing term
if axis == 0
    fx   =  fS.fx ;
    dpdx =  fS.dpdx ;
    
elseif axis == 1
    fx   = fS.fy;
    dpdx =  fS.dpdy ;
    
end

%dpdxApprox = dpdx(x(iB,jB),y(iB,jB),t2);

% dudt + u dudx + v dudy = -dpdx + mu (dudx2 + dudy2)
% if no twilightZone forcing, uxx = 1/mu * ( dp/dx - fx) - uyy
% ignore uyy term for a moment
% all the rest terms are equals to 0 on a no slip wall

if tw == 1
    F = fx(x(iB,jB),y(iB,jB),t2);
else
    F = 0;
end

Uxx = (1./mu).* ( dpdxApprox - F);

if tw > 0
    
    u     = fS.u;
    v     = fS.v;
    if axis == 0
        dudt  = fS.dudt ;
        dudx  = fS.dudx ;
        dudy  = fS.dudy ;
        dudx2 = fS.dudx2 ;
    elseif axis == 1
        dudt  = fS.dvdt;
        dudx  = fS.dvdx;
        dudy  = fS.dvdy;
        dudx2 = fS.dvdy2 ;
    end
    
    Uxx = (1./mu).*(dudt(x(iB,jB),y(iB,jB),t2) ...
        +  u(x(iB,jB),y(iB,jB),t2).*dudx(x(iB,jB),y(iB,jB),t2) ...
        +  v(x(iB,jB),y(iB,jB),t2).*dudy(x(iB,jB),y(iB,jB),t2) ...
        + dpdxApprox ...
        - fx(x(iB,jB),y(iB,jB),t2)) ...
        - dudx2(x(iB,jB),y(iB,jB),t2);
    
elseif tw == 0
    
    Uxx =  Uxx - getUxx(U,iB,jB,axis,h);
    
end

%Uxx = Uxx - (1./mu)*mu4;



end


function dpdx = getDpdx(fS,i,j,count,axis,h,order)

if count == 1
    % at the predictor step
    tSpan = fS.tSpan;
    t2     = tSpan(end);
    t1     = tSpan(end-1) ;
    tp1    = tSpan(end-2) ;
    tp2    = tSpan(end-3) ;
    
    % got this coeffs from miscellaneous.m
    coeff(1)=         -(t2*tp1 + t2*tp2 - tp1*tp2 - t2^2)/...
        ((t1 - tp1)*(t1 - tp2));
    coeff(2)=           (t1*t2 - t1*tp2 + t2*tp2 - t2^2)/...
        ((t1 - tp1)*(tp1 - tp2));
    coeff(3)=-(t1*t2 - t1*tp1 + t2*tp1 - t2^2)/...
        (t1*tp1 - t1*tp2 - tp1*tp2 + tp2^2);
    
    PC = fS.PC;
    PP1 = fS.PP1;
    PP2 = fS.PP2;
    
    PE =  coeff(1)*PC ...
        + coeff(2)*PP1 ...
        + coeff(3)*PP2;
    
elseif count > 1
    
    PE = fS.PN;
    
end

xShift =   (axis==0);
yShift =   (axis==1);

if order == 1
    
    dpdx  = ((-1)* PE(i + xShift*2, j + yShift*2) ...
        + 8* PE(i + xShift, j + yShift) + ...
        (+1)* PE(i - xShift*2, j - yShift*2) ...
        - 8* PE(i - xShift, j - yShift)  ...
        )/(12*h);
    
    
    %     dpdx  = ( ...
    %     + 1* PE(i + xShift, j + yShift) + ...
    %     - 1* PE(i - xShift, j - yShift)  ...
    %     )/(2*h);
    
elseif order == 2
    dpdx  = ((-1)* PE(i + xShift*2, j + yShift*2) ...
        + 16* PE(i + xShift, j + yShift) +   ...
        (-1)* PE(i - xShift*2, j - yShift*2) ...
        + 16* PE(i - xShift, j - yShift) +   ...
        (-30)* PE(i           , j          ) )/(12*h^2);
    
    %     dpdx  = ( ...
    %     + 1* PE(i + xShift, j + yShift) +   ...
    %     + 1* PE(i - xShift, j - yShift) +   ...
    %     (-2)* PE(i           , j          ) )/(h^2);
    
end

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



function approxUx = getCompUx(U,fS,bcx,bcy,axis,side,pos)

tw = fS.tw;
t2 = fS.t2;
x  = fS.x;
y  = fS.y;

% keep the exact ux here for debug
% dudx = (axis==1)*fS.dudx + (axis==0)*fS.dvdy;
% U(bcx,bcy) = dudx(x(iB,jB),y(iB,jB),t2);

% get the index for the corresponding boundary points

iB = (axis==0)*(bcx + (-1)^side*pos) + (axis==1)*bcx;
jB = (axis==1)*(bcy + (-1)^side*pos) + (axis==0)*bcy;

% if twilightZone forcing is given, use exact dudx,
% otherwise, use approximated dudx
if tw > 0
    
    if axis == 0
        dudx  = fS.dvdy;
    elseif axis == 1
        dudx = fS.dudx;
    end
    
    approxUx = -dudx(x(iB,jB),y(iB,jB),t2);
    
elseif tw == 0
    hx = fS.hx;
    hy = fS.hy;
    h    = (axis==0)*hy  + (axis==1)*hx;
    axis = abs(axis-1);
    
    approxUx = -getUx(U,iB,jB,axis,h);
    
end

end


function approxUx = getUx(U,iB,jB,axis,h)

xShift =   (axis==0);

yShift =   (axis==1);

approxUx = ((-1)* U(iB + xShift*2, jB + yShift*2) ...
    + 8* U(iB + xShift, jB + yShift) + ...
    (+1)* U(iB - xShift*2, jB - yShift*2) ...
    - 8* U(iB - xShift, jB - yShift)  ...
    )/(12*h);
end

% function approxUxy = getCompUxy(U,fS,bcx,bcy,axis,side,pos)
%
% tw = fS.tw;
% % keep the exact ux here for debug
% % dudx = (axis==1)*fS.dudx + (axis==0)*fS.dvdy;
% % U(bcx,bcy) = dudx(x(iB,jB),y(iB,jB),t2);
%
% % get the index for the corresponding boundary points
%
% iB = (axis==0)*(bcx + (-1)^side*pos) + (axis==1)*bcx;
% jB = (axis==1)*(bcy + (-1)^side*pos) + (axis==0)*bcy;
%
% % if twilightZone forcing is given, use exact dudx,
% % otherwise, use approximated dudx
% if tw > 0
%     if axis == 0
%         dudx  = fS.dvdy;
%     elseif axis == 1
%         dudx = fS.dudx;
%     end
%
%     approxUx = -dudx(x(iB,jB),y(iB,jB),t2);
%
% elseif tw == 0
%     hx = fS.hx;
%     hy = fS.hy;
%     h    = (axis==0)*hy  + (axis==1)*hx;
%     axis = abs(axis-1);
%
%     approxUx = -getUx(U,iB,jB,axis,h);
%
% end
%
% end

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

function [D1u,D2u,D1v,D2v] = getDs(fS,iB,jB,rx,ry,hx,hy)

tw = fS.tw;
x  = fS.x;
y  = fS.y;
t2  = fS.t2;

if tw == 1
    
    dudx  = fS.dudx;
    dudy  = fS.dudy;
    
    dvdx  = fS.dvdx;
    dvdy  = fS.dvdy;
    
    dudx2 = fS.dudx2;
    dvdx2 = fS.dvdx2;
    
    dudy2 = fS.dudy2;
    dvdy2 = fS.dvdy2;
    
    
    D1u = rx*hx*dudx(x(iB,jB),y(iB,jB),t2) + ry*hy*dudy(x(iB,jB),y(iB,jB),t2);
    
    D2u = (rx*hx)^2/2*dudx2(x(iB,jB),y(iB,jB),t2) + (ry*hy)^2/2*dudy2(x(iB,jB),y(iB,jB),t2) - rx*hx*ry*hy*dvdy2(x(iB,jB),y(iB,jB),t2);
    
    D1v = rx*hx*dvdx(x(iB,jB),y(iB,jB),t2) + ry*hy*dvdy(x(iB,jB),y(iB,jB),t2);
    
    D2v = (rx*hx)^2/2*dvdx2(x(iB,jB),y(iB,jB),t2) + (ry*hy)^2/2*dvdy2(x(iB,jB),y(iB,jB),t2) - rx*hx*ry*hy*dudx2(x(iB,jB),y(iB,jB),t2);
    
    
elseif tw==0
    
    D1u = 0;
    D2u = 0;
    D1v = 0;
    D2v = 0;
    
end

end

function mu4 = getArtifitialDissipation( fS,axis,i,j,hx,hy)

UN = fS.UN;
VN = fS.VN;

ad41 = 2;
ad42 = 2;

graduv1 = 0.25*(abs((-UN(i+2,j) + 8*UN(i+1,j) - 8*UN(i-1,j) + UN(i-2,j))/(12*hx)) + ...
    abs((-UN(i,j+2) + 8*UN(i,j+1) - 8*UN(i,j-1) + UN(i,j-2))/(12*hy))+...
    abs((-VN(i+2,j) + 8*VN(i+1,j) - 8*VN(i-1,j) + VN(i-2,j))/(12*hx))+...
    abs((-VN(i,j+2) + 8*VN(i,j+1) - 8*VN(i,j-1) + VN(i,j-2))/(12*hy)));

maxgrad = max(max(graduv1));

if axis == 1 %axis 1 u is the tangential velocity
    
    lapU = (UN(i+2,j) -4*UN(i+1,j) +6*UN(i,j) -4*UN(i-1,j)+UN(i-2,j)) + ...
        (UN(i,j+2) -4*UN(i,j+1) +6*UN(i,j) -4*UN(i,j-1)+UN(i,j-2));
    
elseif axis ==0
    
    lapU = (VN(i+2,j) -4*VN(i+1,j) +6*VN(i,j) -4*VN(i-1,j)+VN(i-2,j)) + ...
        (VN(i,j+2) -4*VN(i,j+1) +6*VN(i,j) -4*VN(i,j-1)+VN(i,j-2));
        
end

mu4 = (ad41+graduv1*ad42).*(hx^2 + hy^2);
%mu4 = -(ad41+graduv1*ad42).*lapU;

end
