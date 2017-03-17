function [count,U,V,rhsuC,rhsvC,maxgrad] = VelocitySolver(t1,t2,count,fS)
%% read in fS
dt  = fS.dt;
dte = fS.dte;
dtn = fS.dtn;

tOrder    = fS.tOrder;
tMethod   = fS.tMethod;
tExplicit = fS. tExplicit;
imTime    = fS.imTime;
directSolve = fS.directSolve;
twilightZone = fS.twilightZone; 
tw = fS.tw; 

beta = fS.beta;
tref = fS.tref;
g    = fS.g;
mu   = fS.mu;
mubx = fS.mubx;
muby = fS.muby;

x    = fS.x;
y    = fS.y;
Nxg    = fS.Nxg;
Nyg    = fS.Nyg;
hx   = fS.hx;
hy   = fS.hy;

BC   = fS.BC;
ad41 = fS.ad41;
ad42 = fS.ad42;

timeint = fS.preTime;
if count>1
    timeint = fS.corrTime;
end

bcp  = fS.bcp;
ad21 = fS.ad21;
ad22 = fS.ad22;
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
PC = fS.PC;
PP1 = fS.PP1;
PP2 = fS.PP2;
PP3 = fS.PP3;
PP4 = fS.Pextr;

TemN = fS.TemN;

rhsuP1=fS.rhsuP1;
rhsuP2=fS.rhsuP2;
rhsuP3=fS.rhsuP3;

rhsvP1=fS.rhsvP1;
rhsvP2=fS.rhsvP2;
rhsvP3=fS.rhsvP3;

u   =fS.u;
v   =fS.v;
fx  =fS.fx;
fy  =fS.fy;
dvdt=fS.dvdt;
dudt=fS.dudt;
tem =fS.tem;
dvdx=fS.dvdx;
dudy=fS.dudy;
fxE  =fS.fxE;
fyE  =fS.fyE;
fxI  =fS.fxI;
fyI  =fS.fyI;

ia = fS.ia;
ib = fS.ib;
ja = fS.ja;
jb = fS.jb;

%% solve velocity

% initialization
U = zeros(Nxg,Nyg);
V = zeros(Nxg,Nyg);
d4u1 = zeros(Nxg,Nyg);
d4v1 = zeros(Nxg,Nyg);
maxgrad = 0;

% if BC == 2 || BC == 6
%     i = 3:Nxg-3;
%     j = 3:Nyg-3;
% else
%     i = 4:Nxg-3;
%     j = 4:Nyg-3;
% end

i = ia+1:ib-1;
j = ja+1:jb-1;

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

if cons == 1
    if WENO == 1
        dirx = [1,0];
        diry = [0,1];
        
        [UF,VF] = getFaceValues(fS,UN,VN,hx,hy);
        
        ux = WENOFlux(UN,UN,dirx,hx,i,j,UF);
        uy = WENOFlux(UN,VN,diry,hy,i,j,VF);
        
        vx = WENOFlux(VN,UN,dirx,hx,i,j,UF);
        vy = WENOFlux(VN,VN,diry,hy,i,j,VF);
    else
        ux(i,j) = (-UN(i+2,j) + 8*UN(i+1,j) - 8*UN(i-1,j) + UN(i-2,j))/(12*hx);
        uy(i,j) = (-UN(i,j+2) + 8*UN(i,j+1) - 8*UN(i,j-1) + UN(i,j-2))/(12*hy);
        vx(i,j) = (-VN(i+2,j) + 8*VN(i+1,j) - 8*VN(i-1,j) + VN(i-2,j))/(12*hx);
        vy(i,j) = (-VN(i,j+2) + 8*VN(i,j+1) - 8*VN(i,j-1) + VN(i,j-2))/(12*hy);
    end
    
    convu(i,j) = UN(i,j).*ux(i,j) + VN(i,j).*uy(i,j);
    convv(i,j) = UN(i,j).*vx(i,j) + VN(i,j).*vy(i,j);
    
else
    
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

rhsuC =  fxE(x(i,j),y(i,j),t1) ...
    + d4u1(i,j) + ...
    - convu(i,j)...
    - (-PN(i+2,j) + 8*PN(i+1,j) - 8*PN(i-1,j) + PN(i-2,j))/(12*hx);

rhsvC = beta*g*(TemN(i,j)-tref) + fyE(x(i,j),y(i,j),t1) ...
    + d4v1(i,j) + ...
    - convv(i,j) ...
    - (-PN(i,j+2) + 8*PN(i,j+1) - 8*PN(i,j-1) + PN(i,j-2))/(12*hy);


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


if count == 0
    return
end


if tExplicit==1
    
    if tOrder==2
        
        U(i,j) = timeint(1)*rhsuC...
            + timeint(2)*rhsuP1...
            + UC(i,j); %
        
        V(i,j) = timeint(1)*rhsvC...
            + timeint(2)*rhsvP1...
            + VC(i,j);
        
        
    elseif tOrder==4
        if tMethod==1
            U(i,j) = timeint(1)*rhsuC...
                + timeint(2)*rhsuP1...
                + timeint(3)*rhsuP2...
                + timeint(4)*rhsuP3...
                + UC(i,j); %
            
            V(i,j) = timeint(1)*rhsvC...
                + timeint(2)*rhsvP1...
                + timeint(3)*rhsvP2...
                + timeint(4)*rhsvP3...
                + VC(i,j);
            
        elseif tMethod==2
            
            U(i,j) = timeint(1)*rhsuC...
                + timeint(2)*rhsuP1...
                + timeint(3)*rhsuP2...
                + timeint(4)*rhsuP3...
                + 12/25*( 4  *UC(i,j) ...
                -3  *UP1(i,j) ...
                +4/3*UP2(i,j) ...
                -1/4*UP3(i,j)); %
            
            V(i,j) = timeint(1)*rhsvC...
                + timeint(2)*rhsvP1...
                + timeint(3)*rhsvP2...
                + timeint(4)*rhsvP3...
                + 12/25*( 4  *VC(i,j) ...
                -3  *VP1(i,j) ...
                +4/3*VP2(i,j) ...
                -1/4*VP3(i,j));
        end
        
    end
    
    
elseif tExplicit==0
    
    
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
        
        U(i,j) = timeint(1)*rhsuC...
            + timeint(2)*rhsuP1...
            + UC(i,j)...
            + imTime(2)*rhsuImC ...
            + imTime(1)*fxN;
        
        V(i,j) = timeint(1)*rhsvC...
            + timeint(2)*rhsvP1...
            + VC(i,j)...
            + imTime(2)*rhsvImC...
            + imTime(1)*fyN;
                
    elseif tOrder==4
        
        
        if tMethod==1
            
            if tw == 1
                totalFxI = imTime(1)*fxI(x(i,j),y(i,j),t2) ;
                totalFyI = imTime(1)*fyI(x(i,j),y(i,j),t2) ;
                
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
            
            
            U(i,j) = timeint(1)*rhsuC...
                + timeint(2)*rhsuP1...
                + timeint(3)*rhsuP2...
                + timeint(4)*rhsuP3...
                + UC(i,j)...
                + imTime(2)*rhsuImC...
                + imTime(3)*rhsuImP1...
                + imTime(4)*rhsuImP2...
                + totalFxI; %
            
            V(i,j) = timeint(1)*rhsvC...
                + timeint(2)*rhsvP1...
                + timeint(3)*rhsvP2...
                + timeint(4)*rhsvP3...
                + VC(i,j)...
                + imTime(2)*rhsvImC...
                + imTime(3)*rhsvImP1...
                + imTime(4)*rhsvImP2...
                + totalFyI;
            
        elseif tMethod==2
            
            if tw == 1
                
                totalFxI = imTime(1)*fxI(x(i,j),y(i,j),t2);
                totalFyI = imTime(1)*fyI(x(i,j),y(i,j),t2);
            
            elseif tw==0
            
                totalFxI = 0;
                totalFyI = 0;
                
            end
            
            U(i,j) = timeint(1)*rhsuC...
                + timeint(2)*rhsuP1...
                + timeint(3)*rhsuP2...
                + timeint(4)*rhsuP3...
                + 12/25*( 4  *UC(i,j) ...
                -3  *UP1(i,j) ...
                +4/3*UP2(i,j) ...
                -1/4*UP3(i,j))...
                +totalFxI; %
            
            V(i,j) = timeint(1)*rhsvC...
                + timeint(2)*rhsvP1...
                + timeint(3)*rhsvP2...
                + timeint(4)*rhsvP3...
                + 12/25*( 4  *VC(i,j) ...
                -3  *VP1(i,j) ...
                +4/3*VP2(i,j) ...
                -1/4*VP3(i,j))...
                + totalFyI;
        end
        
    end
    
end


%% calculating ghost points
% boundary points
if tExplicit == 1
    switch BC
        case 1
            i = 3;
            U(i,:) = u(x(i,1),y(i,:),t2);
            V(i,:) = v(x(i,1),y(i,:),t2);
            
            i = Nxg-2;
            U(i,:) = u(x(i,1),y(i,:),t2);
            V(i,:) = v(x(i,1),y(i,:),t2);
            
            j = 3;
            U(:,j) = u(x(:,j),y(1,j),t2);
            V(:,j) = v(x(:,j),y(1,j),t2);
            
            j = Nyg-2;
            U(:,j) = u(x(:,j),y(1,j),t2);
            V(:,j) = v(x(:,j),y(1,j),t2);
            %x-ghost points
            
            for k = 1:2
                %ghost v
                i = bcp(1,k); %3; Nxg-2
                j = 4:Nyg - 3;
                index = (-1)^(k-1);
                Vx(:,j) = V(i+index*(1:4),j);
                Ux(:,j) = U(i+index*(1:4),j);
                
                gradbcx = 1/4*( ...
                    2*abs((-U(i,j+2) + 8*U(i,j+1) - 8*U(i,j-1) + U(i,j-2))/(12*hy))+...
                    abs((-V(i,j+2) + 8*V(i,j+1) - 8*V(i,j-1) + V(i,j-2))/(12*hy))+ ...
                    abs((-3*V(i,j) + 4*Vx(1,j) - Vx(2,j))/(2*hx)));
                
                lapvbc = (V(i,j+2) -4*V(i,j+1) +6*V(i,j) -4*V(i,j-1)+V(i,j-2));
                d4vbc(i,j) = -(ad41+gradbcx*ad42);
                d2vbc(i,j) =  (ad21+gradbcx*ad22);
                
                patch11 = mu/hx^2 * eye(length(4:Nyg - 3)) + diag(d4vbc(i,j))*(-4) + diag(d2vbc(i,j));
                %patch11 = eye(length(4:Nyg - 3));    % !!!!!
                patch12 = diag(d4vbc(i,j));
                patch21 = -6*eye(length(4:Nyg - 3));
                patch22 =    eye(length(4:Nyg - 3));
                
                AVbc = [patch11 patch12;patch21 patch22];
                
                if count == 1
                    
                    %if tOrder==2
                    dpdyt = + 3*(PC(i,j+1) - PC(i,j-1))/(2*hy) ...
                        - 3*(PP1(i,j+1) - PP1(i,j-1))/(2*hy) ...
                        +   (PP2(i,j+1) - PP2(i,j-1))/(2*hy);
                    %
                    dpdya = ((((PC(i,j+1) - PC(i,j-1))/(2*hy)  ...
                        - (PP1(i,j+1) - PP1(i,j-1))/(2*hy))/dt ...
                        - ((PP1(i,j+1) - PP1(i,j-1))/(2*hy) ...
                        - (PP2(i,j+1) - PP2(i,j-1))/(2*hy))/dte) ...
                        * (dtn+dt)/(dt+dte) ...
                        + ((PC(i,j+1) - PC(i,j-1))/(2*hy)  ...
                        - (PP1(i,j+1) - PP1(i,j-1))/(2*hy))/dt)*dtn ...
                        + (PC(i,j+1) - PC(i,j-1))/(2*hy);
                    %                 elseif tOrder==4
                    %
                    %                     dpdya = + 5*(PC(i,j+1) - PC(i,j-1))/(2*hy) ...
                    %                         - 10*(PP1(i,j+1) - PP1(i,j-1))/(2*hy) ...
                    %                         + 10*(PP2(i,j+1) - PP2(i,j-1))/(2*hy) ...
                    %                         - 5*(PP3(i,j+1) - PP3(i,j-1))/(2*hy) ...
                    %                         + 1*(PP4(i,j+1) - PP4(i,j-1))/(2*hy);
                    
                    %end
                    
                    %dpdya = dpdy(x(i,j),y(i,j),t2);
                else
                    dpdya = (PN(i,j+1) - PN(i,j-1))/(2*hy);
                    %dpdya = dpdy(x(i,j),y(i,j),t2);
                end
                
                bVbc1 = dvdt(x(i,j),y(i,j),t2) - beta*g*(tem(x(i,j),y(i,j),t2)-tref) ...
                    - fy(x(i,j),y(i,j),t2) ... - d4vbc(i,j).*(Vx(2,j) -4*Vx(1,j) +6*V(i,j)) ...
                    - d4vbc(i,j).*lapvbc ...
                    + u(x(i,j),y(i,j),t2).*dvdx(x(i,j),y(i,j),t2)...
                    + V(i,j).*(V(i,j+1) - V(i,j-1))/(2*hy)...
                    - (mu).*(V(i,j+1) -2*V(i,j) + V(i,j-1))/(hy^2)...
                    - d2vbc(i,j).*(V(i,j+1) -2*V(i,j) + V(i,j-1))...
                    + dpdya ...
                    - (mu).*(Vx(1,j)-2*V(i,j))/hx^2 ...
                    - d2vbc(i,j).*(Vx(1,j)-2*V(i,j));
                
                %bVbc1 = -(Vx(1,j)-2*V(i,j)); %!!!!!
                
                
                bVbc1 = bVbc1 - d4vbc(i,j).*(Vx(2,j) - 4*Vx(1,j) + 6*V(i,j));
                
                
                bVbc2 = -(+15*V(i,j)-20*Vx(1,j) + 15*Vx(2,j) -6*Vx(3,j) + Vx(4,j));
                bVbc = [bVbc1 bVbc2];
                
                vtemp = AVbc\bVbc';
                
                V(i-index*(1),j) = vtemp(1:length(4:Nyg - 3));
                V(i-index*(2),j) = vtemp(length(4:Nyg - 3)+1:end);
                
                V(i-index*(1),2) = v(x(i-index*(1),2),y(i-index*(1),2),t2);
                V(i-index*(2),2) = v(x(i-index*(2),2),y(i-index*(2),2),t2);
                V(i-index*(1),Nyg-1) = v(x(i-index*(1),Nyg-1),y(i-index*(1),Nyg-1),t2);
                V(i-index*(2),Nyg-1) = v(x(i-index*(2),Nyg-1),y(i-index*(2),Nyg-1),t2);
                %ghost u
                A = [16 -1;index*(-8) index*(1)];
                
                for  j = 5:Nyg - 4;
                    
                    rhs = -(-V(i+2,j+2) + 8*V(i+1,j+2) - 8*V(i-1,j+2) + V(i-2,j+2)) ...
                        +8*(-V(i+2,j+1) + 8*V(i+1,j+1) - 8*V(i-1,j+1) + V(i-2,j+1)) ...
                        -8*(-V(i+2,j-1) + 8*V(i+1,j-1) - 8*V(i-1,j-1) + V(i-2,j-1)) ...
                        + (-V(i+2,j-2) + 8*V(i+1,j-2) - 8*V(i-1,j-2) + V(i-2,j-2));
                    
                    b(1) = -rhs/(12*hx*12*hy)*12*hx^2 - (-Ux(2,j) + 16*Ux(1,j) -30*U(i,j));
                    
                    b(2) = -(-V(i,j+2) + 8*V(i,j+1) - 8*V(i,j-1) + V(i,j-2))/(12*hy)*12*hx...
                        - index*(-Ux(2,j) + 8*Ux(1,j));
                    x1 = A\b';
                    U(i-index*(1),j) = x1(1);
                    U(i-index*(2),j) = x1(2);
                end
            end
            
            % ghost points in y direction
            for k = 1:2
                %ghost v
                j = bcp(2,k);%3; Nyg-2
                i = 4:Nxg-3;
                index = (-1)^(k-1);
                Vy(i,:) = V(i,j+index*(1:4));
                Uy(i,:) = U(i,j+index*(1:4));
                
                gradbcy = 0.25*(2*abs((-U(i+2,j) + 8*U(i+1,j) - 8*U(i-1,j) + U(i-2,j))/(12*hx)) + ...
                    abs((-3*U(i,j) + 4*Uy(i,1) - Uy(i,2))/(2*hy))+...
                    abs((-V(i+2,j) + 8*V(i+1,j) - 8*V(i-1,j) + V(i-2,j))/(12*hx)));
                
                lapubc = U(i+2,j) -4*U(i+1,j) +6*U(i,j) -4*U(i-1,j)+U(i-2,j);
                
                d4ubc(i,j) = -(ad41+gradbcy*ad42);
                d2ubc(i,j) =  (ad21+gradbcy*ad22);
                
                
                patch11 = mu/hy^2 * eye(length(4:Nyg - 3)) + diag(d4ubc(i,j))*(-4) + diag(d2ubc(i,j));
                %patch11 = eye(length(4:Nyg - 3)) ;%!!!!!!
                patch12 = diag(d4ubc(i,j));
                patch21 = -6*eye(length(4:Nyg - 3));
                patch22 =    eye(length(4:Nyg - 3));
                
                AUbc = [patch11 patch12;patch21 patch22];
                
                if count == 1
                    % if tOrder==2
                    %                 dpdx = - 3*(-(P1(i+1,j) - P1(i-1,j)))/(2*hx) ...
                    %                     + 3*(-(P0(i+1,j) - P0(i-1,j)))/(2*hx) ...
                    %                     -   (-(Pextr(i+1,j) - Pextr(i-1,j)))/(2*hx);
                    %
                    dpdxa = ((((PC(i+1,j) - PC(i-1,j))/(2*hx)  ...
                        - (PP1(i+1,j) - PP1(i-1,j))/(2*hx))/dt ...
                        - ((PP1(i+1,j) - PP1(i-1,j))/(2*hx) ...
                        - (PP2(i+1,j) - PP2(i-1,j))/(2*hx))/dte) ...
                        * (dtn+dt)/(dt+dte) ...
                        + ((PC(i+1,j) - PC(i-1,j))/(2*hx)  ...
                        - (PP1(i+1,j) - PP1(i-1,j))/(2*hx))/dt)*dtn ...
                        + (PC(i+1,j) - PC(i-1,j))/(2*hx);
                    
                    %dpdxa = dpdx(x(i,j),y(i,j),t2);
                    
                    %                 elseif tOrder==4
                    %
                    %                     dpdxa = + 5*(PC(i+1,j) - PC(i-1,j))/(2*hx) ...
                    %                         - 10*(PP1(i+1,j) - PP1(i-1,j))/(2*hx) ...
                    %                         + 10*(PP2(i+1,j) - PP2(i-1,j))/(2*hx) ...
                    %                         - 5*(PP3(i+1,j) - PP3(i-1,j))/(2*hx) ...
                    %                         + 1*(PP4(i+1,j) - PP4(i-1,j))/(2*hx);
                    
                    % end
                    
                else
                    dpdxa = (PN(i+1,j) - PN(i-1,j))/(2*hx);
                    %                  dpdxa = dpdx(x(i,j),y(i,j),t2);
                end
                
                
                bUbc1 = dudt(x(i,j),y(i,j),t2) ...
                    - fx(x(i,j),y(i,j),t2) ...- d4ubc(i,j).*(Uy(i,2) -4*Uy(i,1) +6*U(i,j)) ...
                    - d4ubc(i,j).*lapubc ...
                    + v(x(i,j),y(i,j),t2).*dudy(x(i,j),y(i,j),t2) ...
                    + U(i,j).*(U(i+1,j) - U(i-1,j))/(2*hx) ...
                    + dpdxa ...
                    - mu*(U(i+1,j) -2*U(i,j) + U(i-1,j))/(hx^2)...
                    - d2ubc(i,j).*(U(i+1,j) -2*U(i,j) + U(i-1,j))...
                    - mu*(Uy(i,1)-2*U(i,j))/hy^2 ...
                    - d2ubc(i,j).*(Uy(i,1)-2*U(i,j));
                %bUbc1 = - (Uy(i,1)-2*U(i,j)); %!!!!
                
                
                bUbc1 = bUbc1 - d4ubc(i,j).*(Uy(i,2) -4*Uy(i,1) +6*U(i,j));
                
                bUbc2 = -(+15*U(i,j)-20*Uy(i,1) + 15*Uy(i,2) -6*Uy(i,3) + Uy(i,4));
                bUbc = [bUbc1' bUbc2'];
                
                utemp = AUbc\bUbc';
                U(i,j-index*(1)) = utemp(1:length(4:Nyg - 3));
                U(i,j-index*(2)) = utemp(length(4:Nyg - 3)+1:end);
                
                
                % fill in corners that needed in the corss derivative
                U(2,j-index*(1)) = u(x(2,j-index*(1)),y(2,j-index*(1)),t2);
                U(2,j-index*(2)) = u(x(2,j-index*(2)),y(2,j-index*(2)),t2);
                U(Nxg-1,j-index*(1)) = u(x(Nxg-1,j-index*(1)),y(Nxg-1,j-index*(1)),t2);
                U(Nxg-1,j-index*(2)) = u(x(Nxg-1,j-index*(2)),y(Nxg-1,j-index*(2)),t2);
                
                A = [16 -1;index*(-8) index*1];
                for  i = 5:Nxg - 4;
                    
                    rhs = -(-U(i+2,j+2) + 8*U(i+1,j+2) - 8*U(i-1,j+2) + U(i-2,j+2)) ...
                        +8*(-U(i+2,j+1) + 8*U(i+1,j+1) - 8*U(i-1,j+1) + U(i-2,j+1)) ...
                        -8*(-U(i+2,j-1) + 8*U(i+1,j-1) - 8*U(i-1,j-1) + U(i-2,j-1)) ...
                        + (-U(i+2,j-2) + 8*U(i+1,j-2) - 8*U(i-1,j-2) + U(i-2,j-2));
                    
                    b(1) = -rhs/(12*hx*12*hy)*12*hy^2 - (-Vy(i,2) + 16*Vy(i,1) -30*V(i,j));
                    b(2) = -(-U(i+2,j) + 8*U(i+1,j) - 8*U(i-1,j) + U(i-2,j))/(12*hx)*12*hy...
                        - index*(-Vy(i,2) + 8*Vy(i,1));
                    x1 = A\b';
                    V(i,j-index*(1)) = x1(1);
                    V(i,j-index*(2)) = x1(2);
                end
                
            end
            %% Calculate 'corner' ghost points
            %lower left corner
            A = [16 -1 8/(12*hx*12*hy)*12*hx^2 0;-8 1 0 0; ...
                8/(12*hx*12*hy)*12*hy^2 0 16 -1;0 0 -8 1];
            
            i =3;
            j =4;
            rhs = -(-V(i+2,j+2) + 8*V(i+1,j+2) - 8*V(i-1,j+2) + V(i-2,j+2)) ...
                +8*(-V(i+2,j+1) + 8*V(i+1,j+1) - 8*V(i-1,j+1) + V(i-2,j+1)) ...
                -8*(-V(i+2,j-1) + 8*V(i+1,j-1) - 8*V(i-1,j-1) + V(i-2,j-1)) ...
                + (-V(i+2,j-2) - 8*V(i-1,j-2) + V(i-2,j-2));
            
            bc(1) = -rhs/(12*hx*12*hy)*12*hx^2 - (-U(i+2,j) + 16*U(i+1,j) -30*U(i,j));
            bc(2) = -(-V(i,j+2) + 8*V(i,j+1) - 8*V(i,j-1) + V(i,j-2))/(12*hy)*12*hx...
                - (-U(i+2,j) + 8*U(i+1,j));
            
            j = 3;
            i = 4;
            rhs = -(-U(i+2,j+2) + 8*U(i+1,j+2) - 8*U(i-1,j+2) + U(i-2,j+2)) ...
                +8*(-U(i+2,j+1) + 8*U(i+1,j+1) - 8*U(i-1,j+1) ) ...
                -8*(-U(i+2,j-1) + 8*U(i+1,j-1) - 8*U(i-1,j-1) + U(i-2,j-1)) ...
                + (-U(i+2,j-2) + 8*U(i+1,j-2) - 8*U(i-1,j-2) + U(i-2,j-2));
            
            bc(3) = -rhs/(12*hx*12*hy)*12*hy^2 - (-V(i,j+2) + 16*V(i,j+1) -30*V(i,j));
            bc(4) = -(-U(i+2,j) + 8*U(i+1,j) - 8*U(i-1,j) + U(i-2,j))/(12*hx)*12*hy...
                - (-V(i,j+2) + 8*V(i,j+1));
            
            x1 = A\bc';
            
            V(i,j-1) = x1(3);
            V(i,j-2) = x1(4);
            i =3;
            j =4;
            U(i-1,j) = x1(1);
            U(i-2,j) = x1(2);
            
            %upper left corner
            
            A = [16 -1 -8/(12*hx*12*hy)*12*hx^2 0;-8 1 0 0; ...
                -8/(12*hx*12*hy)*12*hy^2 0 16 -1;0 0 8 -1];
            
            %         A = [16 -1 0 0;-8 1 0 0; ...
            %             0 0 16 -1;0 0 8 -1];
            
            i = 3;
            j = Nyg-3;
            %         V2(i+1,j+2) = v(x(i+1,j+2),y(i+1,j+2),t2);
            rhs = -(-V(i+2,j+2)  - 8*V(i-1,j+2) + V(i-2,j+2)) ...
                +8*(-V(i+2,j+1) + 8*V(i+1,j+1) - 8*V(i-1,j+1) + V(i-2,j+1)) ...
                -8*(-V(i+2,j-1) + 8*V(i+1,j-1) - 8*V(i-1,j-1) + V(i-2,j-1)) ...
                + (-V(i+2,j-2) +8*V(i+1,j-2)- 8*V(i-1,j-2) + V(i-2,j-2));
            
            bc(1) = -rhs/(12*hx*12*hy)*12*hx^2 - (-U(i+2,j) + 16*U(i+1,j) -30*U(i,j));
            bc(2) = -(-V(i,j+2) + 8*V(i,j+1) - 8*V(i,j-1) + V(i,j-2))/(12*hy)*12*hx...
                - (-U(i+2,j) + 8*U(i+1,j));
            
            j = Nyg-2;
            i = 4;
            %         U2(i-2,j-1) = u(x(i-2,j-1),y(i-2,j-1),t2);
            rhs = -(-U(i+2,j+2) + 8*U(i+1,j+2) - 8*U(i-1,j+2) + U(i-2,j+2)) ...
                +8*(-U(i+2,j+1) + 8*U(i+1,j+1) - 8*U(i-1,j+1) + U(i-2,j+1)) ...
                -8*(-U(i+2,j-1) + 8*U(i+1,j-1) - 8*U(i-1,j-1) ) ...
                + (-U(i+2,j-2) + 8*U(i+1,j-2) - 8*U(i-1,j-2) + U(i-2,j-2));
            
            bc(3) = -rhs/(12*hx*12*hy)*12*hy^2 - (-V(i,j-2) + 16*V(i,j-1) -30*V(i,j));
            bc(4) = -(-U(i+2,j) + 8*U(i+1,j) - 8*U(i-1,j) + U(i-2,j))/(12*hx)*12*hy...
                - (-8*V(i,j-1) + V(i,j-2));
            
            x1 = A\bc';
            
            V(i,j+1) = x1(3);
            V(i,j+2) = x1(4);
            
            i =3;
            j =Nyg-3;
            U(i-1,j) = x1(1);
            U(i-2,j) = x1(2);
            
            
            %uper right corner
            
            A = [16 -1 8/(12*hx*12*hy)*12*hx^2 0;8 -1 0 0; ...
                8/(12*hx*12*hy)*12*hy^2 0 16 -1;0 0 8 -1];
            
            %         A = [16 -1 0 0;-8 1 0 0; ...
            %             0 0 16 -1;0 0 8 -1];
            
            i = Nxg-2;
            j = Nyg-3;
            %         V2(i+1,j+2) = v(x(i+1,j+2),y(i+1,j+2),t2);
            rhs = -(-V(i+2,j+2) + 8*V(i+1,j+2)  + V(i-2,j+2)) ...
                +8*(-V(i+2,j+1) + 8*V(i+1,j+1) - 8*V(i-1,j+1) + V(i-2,j+1)) ...
                -8*(-V(i+2,j-1) + 8*V(i+1,j-1) - 8*V(i-1,j-1) + V(i-2,j-1)) ...
                + (-V(i+2,j-2) +8*V(i+1,j-2)- 8*V(i-1,j-2) + V(i-2,j-2));
            
            bc(1) = -rhs/(12*hx*12*hy)*12*hx^2 - (-U(i-2,j) + 16*U(i-1,j) -30*U(i,j));
            bc(2) = -(-V(i,j+2) + 8*V(i,j+1) - 8*V(i,j-1) + V(i,j-2))/(12*hy)*12*hx...
                - (-8*U(i-1,j) + U(i-2,j));
            
            j = Nyg-2;
            i = Nxg-3;
            %         U2(i-2,j-1) = u(x(i-2,j-1),y(i-2,j-1),t2);
            rhs = -(-U(i+2,j+2) + 8*U(i+1,j+2) - 8*U(i-1,j+2) + U(i-2,j+2)) ...
                +8*(-U(i+2,j+1) + 8*U(i+1,j+1) - 8*U(i-1,j+1) + U(i-2,j+1)) ...
                -8*( + 8*U(i+1,j-1) - 8*U(i-1,j-1) + U(i-2,j-1)) ...
                + (-U(i+2,j-2) + 8*U(i+1,j-2) - 8*U(i-1,j-2) + U(i-2,j-2));
            
            bc(3) = -rhs/(12*hx*12*hy)*12*hy^2 - (-V(i,j-2) + 16*V(i,j-1) -30*V(i,j));
            bc(4) = -(-U(i+2,j) + 8*U(i+1,j) - 8*U(i-1,j) + U(i-2,j))/(12*hx)*12*hy...
                - (-8*V(i,j-1) + V(i,j-2));
            
            x1 = A\bc';
            
            V(i,j+1) = x1(3);
            V(i,j+2) = x1(4);
            
            i = Nxg-2;
            j =Nyg-3;
            U(i+1,j) = x1(1);
            U(i+2,j) = x1(2);
            
            
            %lower right corner
            
            A = [16 -1 -8/(12*hx*12*hy)*12*hx^2 0;8 -1 0 0; ...
                -8/(12*hx*12*hy)*12*hy^2 0 16 -1;0 0 -8 1];
            
            %         A = [16 -1 0 0;-8 1 0 0; ...
            %             0 0 16 -1;0 0 8 -1];
            
            i = Nxg-2;
            j = 4;
            %         V2(i+1,j+2) = v(x(i+1,j+2),y(i+1,j+2),t2);
            rhs = -(-V(i+2,j+2) + 8*V(i+1,j+2) - 8*V(i-1,j+2) + V(i-2,j+2)) ...
                +8*(-V(i+2,j+1) + 8*V(i+1,j+1) - 8*V(i-1,j+1) + V(i-2,j+1)) ...
                -8*(-V(i+2,j-1) + 8*V(i+1,j-1) - 8*V(i-1,j-1) + V(i-2,j-1)) ...
                + (-V(i+2,j-2) +8*V(i+1,j-2) + V(i-2,j-2));
            
            bc(1) = -rhs/(12*hx*12*hy)*12*hx^2 - (-U(i-2,j) + 16*U(i-1,j) -30*U(i,j));
            bc(2) = -(-V(i,j+2) + 8*V(i,j+1) - 8*V(i,j-1) + V(i,j-2))/(12*hy)*12*hx...
                - (-8*U(i-1,j) + U(i-2,j));
            
            j = 3;
            i = Nxg-3;
            %         U2(i-2,j-1) = u(x(i-2,j-1),y(i-2,j-1),t2);
            rhs = -(-U(i+2,j+2) + 8*U(i+1,j+2) - 8*U(i-1,j+2) + U(i-2,j+2)) ...
                +8*( + 8*U(i+1,j+1) - 8*U(i-1,j+1) + U(i-2,j+1)) ...
                -8*(-U(i+2,j-1) + 8*U(i+1,j-1) - 8*U(i-1,j-1) + U(i-2,j-1)) ...
                + (-U(i+2,j-2) + 8*U(i+1,j-2) - 8*U(i-1,j-2) + U(i-2,j-2));
            
            bc(3) = -rhs/(12*hx*12*hy)*12*hy^2 - (-V(i,j+2) + 16*V(i,j+1) -30*V(i,j));
            bc(4) = -(-U(i+2,j) + 8*U(i+1,j) - 8*U(i-1,j) + U(i-2,j))/(12*hx)*12*hy...
                - (-V(i,j+2) + 8*V(i,j+1));
            
            x1 = A\bc';
            
            V(i,j-1) = x1(3);
            V(i,j-2) = x1(4);
            
            i = Nxg-2;
            j =4;
            U(i+1,j) = x1(1);
            U(i+2,j) = x1(2);
            
        case 2
            
            V(bcp(1,2),:) = V(bcp(1,1),:);
            U(bcp(1,2),:) = U(bcp(1,1),:);
            V(:,bcp(2,2)) = V(:,bcp(2,1));
            U(:,bcp(2,2)) = U(:,bcp(2,1));
            
            V(bcp(1,1)-1,:) = V(bcp(1,2)-1,:);
            V(bcp(1,1)-2,:) = V(bcp(1,2)-2,:);
            V(bcp(1,2)+1,:) = V(bcp(1,1)+1,:);
            V(bcp(1,2)+2,:) = V(bcp(1,1)+2,:);
            
            U(bcp(1,1)-1,:) = U(bcp(1,2)-1,:);
            U(bcp(1,1)-2,:) = U(bcp(1,2)-2,:);
            U(bcp(1,2)+1,:) = U(bcp(1,1)+1,:);
            U(bcp(1,2)+2,:) = U(bcp(1,1)+2,:);
            
            V(:,bcp(2,1)-1) = V(:,bcp(2,2)-1);
            V(:,bcp(2,1)-2) = V(:,bcp(2,2)-2);
            V(:,bcp(2,2)+1) = V(:,bcp(2,1)+1);
            V(:,bcp(2,2)+2) = V(:,bcp(2,1)+2);
            
            U(:,bcp(2,1)-1) = U(:,bcp(2,2)-1);
            U(:,bcp(2,1)-2) = U(:,bcp(2,2)-2);
            U(:,bcp(2,2)+1) = U(:,bcp(2,1)+1);
            U(:,bcp(2,2)+2) = U(:,bcp(2,1)+2);
            
        case 3
            
            i = 1:3;
            U(i,:) = u(x(i,:),y(i,:),t2);
            V(i,:) = v(x(i,:),y(i,:),t2);
            
            i = Nxg-2:Nxg;
            U(i,:) = u(x(i,:),y(i,:),t2);
            V(i,:) = v(x(i,:),y(i,:),t2);
            
            j = 1:3;
            U(:,j) = u(x(:,j),y(:,j),t2);
            V(:,j) = v(x(:,j),y(:,j),t2);
            
            j = Nyg-2:Nyg;
            U(:,j) = u(x(:,j),y(:,j),t2);
            V(:,j) = v(x(:,j),y(:,j),t2);
            
            
    end
    
    %fill the other corner points
    i = 1:2;
    j = 1:2;
    U(i,j) =u(x(i,j),y(i,j),t2);
    V(i,j) =v(x(i,j),y(i,j),t2);
    j = Nyg-1:Nyg;
    U(i,j) =u(x(i,j),y(i,j),t2);
    V(i,j) =v(x(i,j),y(i,j),t2);
    i = Nxg-1:Nxg;
    j = 1:2;
    U(i,j) =u(x(i,j),y(i,j),t2);
    V(i,j) =v(x(i,j),y(i,j),t2);
    j = Nyg-1:Nyg;
    U(i,j) =u(x(i,j),y(i,j),t2);
    V(i,j) =v(x(i,j),y(i,j),t2);
    
elseif tExplicit==0
    if BC == 4
        J  = 1:Nyg;
        for i = 1:3
            U(i,J) =u(x(i,J),y(i,J),t2);
            V(i,J) =v(x(i,J),y(i,J),t2);
        end
    
        for i = Nxg-2:Nxg
            U(i,J) =u(x(i,J),y(i,J),t2);
            V(i,J) =v(x(i,J),y(i,J),t2);
        end
    
        I  = 1:Nxg;
        for j = 1:3
            U(I,j) =u(x(I,j),y(I,j),t2);
            V(I,j) =v(x(I,j),y(I,j),t2);
        end
        
        for j = Nyg-2:Nyg
            U(I,j) =u(x(I,j),y(I,j),t2);
            V(I,j) =v(x(I,j),y(I,j),t2);
        end
        
        dir=0;
        U = solveIm(Nxg,Nyg,U,dir,directSolve);
        
        dir=1;
        V = solveIm(Nxg,Nyg,V,dir,directSolve);
        
    elseif BC == 2
         
        dir=0;
        U = solveIm(Nxg,Nyg,U,dir,directSolve);
        
        dir=1;
        V = solveIm(Nxg,Nyg,V,dir,directSolve);
        
    elseif BC == 1 || BC == 6
        
        if fS.uvSwitch == 1  % switch solve u first then correct with v
            preIndex     = 1;
            correctIndex = 1;
        elseif fS.uvSwitch == 0 % dont switch solve u first 
            preIndex     = 100;
            correctIndex = 100;
            
        elseif fS.uvSwitch == 2 % dont switch solve v first
            preIndex     = 0;
            correctIndex = 0;
            
        end
        
        %if BC == 1
        if count <= preIndex
            U = solveGhostU(fS,BC,Nxg,Nyg,t2,u,v,x,y,hx,hy,dt,dte,dtn,...
                U,...
                VN,VC,VP1,VP2,...
                PN,PC,PP1,PP2,...
                count,tOrder);
            
            if BC == 1
                
                dir=0;
                U = solveIm(Nxg,Nyg,U,dir,directSolve);
                
            elseif BC == 6
                
                i = 1;
                J = 1:Nyg;
                U(i,J) = 0;
                
                i = 2;
                U(i,J) = 0;
                
                for i = Nxg-2:Nxg
                    U(i,J) = 0;
                end
                
                dir=0;
                U = solveIm(Nxg,Nyg,U,dir,directSolve);
            end
            
            
            V = solveGhostV(V,U,UC,UP1,UP2);
            
            if BC == 1
                
                dir=1;
                V = solveIm(Nxg,Nyg,V,dir,directSolve);
                
            elseif BC == 6
                
                i = 1;
                J = 1:Nyg;
                V(i,J) = 0;
                
                i = 2;
                V(i,J) = 0;
                
                for i = Nxg-2:Nxg
                    V(i,J) = 0;
                end
                
                dir=1;
                V = solveIm(Nxg,Nyg,V,dir,directSolve);
                
                
            end
            
        elseif count > correctIndex
            
            V = solveGhostV(V,UN,UC,UP1,UP2);
            
            if BC == 1
                
                dir=1;
                V = solveIm(Nxg,Nyg,V,dir,directSolve);
                
            elseif BC == 6
                
                i = 1;
                J = 1:Nyg;
                V(i,J) = 0;
                
                i = 2;
                V(i,J) = 0;
                
                for i = Nxg-2:Nxg
                    V(i,J) = 0;
                end
                
                dir=1;
                V = solveIm(Nxg,Nyg,V,dir,directSolve);
                
            end
            
            U = solveGhostU(fS,BC,Nxg,Nyg,t2,u,v,x,y,hx,hy,dt,dte,dtn,...
                U,...
                V,VC,VP1,VP2,...
                PN,PC,PP1,PP2,...
                count,tOrder);
            
            if BC == 1
                
                dir=0;
                U = solveIm(Nxg,Nyg,U,dir,directSolve);
                
            elseif BC == 6
                
                i = 1;
                J = 1:Nyg;
                U(i,J) = 0;
                
                i = 2;
                U(i,J) = 0;
                
                for i = Nxg-2:Nxg
                    U(i,J) = 0;
                end
                
                dir=0;
                U = solveIm(Nxg,Nyg,U,dir,directSolve);
            end
        end
        
    end
    
end

count = count + 1;


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

    function F = WENOFlux(f,fd,dir,h,i,j,Ud)
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

    function W = solveIm(Nxg,Nyg,U,dir,directSolve)
        %% reshape to the rhs vector
        rhsu =zeros(Nxg*Nyg,1);
        
        for l = 1:Nxg
            
            index = ((l-1)*Nyg+1:l*Nyg);
            J     = 1:Nyg;
            rhsu(index) = U(l,J);
            
        end
        
        %% get the LHS matrix
        if directSolve == 0
            
            if dir == 0
                Ll = fS.Lul;
                Lu = fS.Luu;
                p  = fS.Up;
            elseif dir ==1
                Ll = fS.Lvl;
                Lu = fS.Lvu;
                p  = fS.Vp;
            end
            yt = Ll\(p*rhsu);
            Ut  = Lu\yt;
            
        elseif directSolve == 1
            
            if dir == 0
                
                Lu = fS.Lu;
                Ut  = Lu\rhsu;
                
            elseif dir ==1
                
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

    function dp = approximateDpdx(i,j,count,PN,PC,PP1,PP2,h,dt,dte,dtn,dir)
        
        if dir == 0
            if count == 1
                % if tOrder==2
                %                 dpdx = - 3*(-(P1(i+1,j) - P1(i-1,j)))/(2*hx) ...
                %                     + 3*(-(P0(i+1,j) - P0(i-1,j)))/(2*hx) ...
                %                     -   (-(Pextr(i+1,j) - Pextr(i-1,j)))/(2*hx);
                %
                dp = ((((PC(i+1,j) - PC(i-1,j))/(2*h)  ...
                    - (PP1(i+1,j) - PP1(i-1,j))/(2*h))/dt ...
                    - ((PP1(i+1,j) - PP1(i-1,j))/(2*h) ...
                    - (PP2(i+1,j) - PP2(i-1,j))/(2*h))/dte) ...
                    * (dtn+dt)/(dt+dte) ...
                    + ((PC(i+1,j) - PC(i-1,j))/(2*h)  ...
                    - (PP1(i+1,j) - PP1(i-1,j))/(2*h))/dt)*dtn ...
                    + (PC(i+1,j) - PC(i-1,j))/(2*h);
                
                %dpdxa = dpdx(x(i,j),y(i,j),t2);
                
                %                 elseif tOrder==4
                %
                %                     dpdxa = + 5*(PC(i+1,j) - PC(i-1,j))/(2*hx) ...
                %                         - 10*(PP1(i+1,j) - PP1(i-1,j))/(2*hx) ...
                %                         + 10*(PP2(i+1,j) - PP2(i-1,j))/(2*hx) ...
                %                         - 5*(PP3(i+1,j) - PP3(i-1,j))/(2*hx) ...
                %                         + 1*(PP4(i+1,j) - PP4(i-1,j))/(2*hx);
                
                % end
                
            else
                dp = (PN(i+1,j) - PN(i-1,j))/(2*h);
                %                  dpdxa = dpdx(x(i,j),y(i,j),t2);
            end
            
        elseif dir ==1
            if count == 1
                % if tOrder==2
                %                 dpdx = - 3*(-(P1(i+1,j) - P1(i-1,j)))/(2*hx) ...
                %                     + 3*(-(P0(i+1,j) - P0(i-1,j)))/(2*hx) ...
                %                     -   (-(Pextr(i+1,j) - Pextr(i-1,j)))/(2*hx);
                %
                dp = ((((PC(i,j+1) - PC(i,j-1))/(2*h)  ...
                    - (PP1(i,j+1) - PP1(i,j-1))/(2*h))/dt ...
                    - ((PP1(i,j+1) - PP1(i,j-1))/(2*h) ...
                    - (PP2(i,j+1) - PP2(i,j-1))/(2*h))/dte) ...
                    * (dtn+dt)/(dt+dte) ...
                    + ((PC(i,j+1) - PC(i,j-1))/(2*h)  ...
                    - (PP1(i,j+1) - PP1(i,j-1))/(2*h))/dt)*dtn ...
                    + (PC(i,j+1) - PC(i,j-1))/(2*h);
                
                %dpdxa = dpdx(x(i,j),y(i,j),t2);
                
                %                 elseif tOrder==4
                %
                %                     dpdxa = + 5*(PC(i+1,j) - PC(i-1,j))/(2*hx) ...
                %                         - 10*(PP1(i+1,j) - PP1(i-1,j))/(2*hx) ...
                %                         + 10*(PP2(i+1,j) - PP2(i-1,j))/(2*hx) ...
                %                         - 5*(PP3(i+1,j) - PP3(i-1,j))/(2*hx) ...
                %                         + 1*(PP4(i+1,j) - PP4(i-1,j))/(2*hx);
                
                % end
                
            else
                dp = (PN(i,j+1) - PN(i,j-1))/(2*h);
                %                  dpdxa = dpdx(x(i,j),y(i,j),t2);
            end
            
            
        end
    end

    function dudxy = approximateDudxy(i,j,U,hx,hy)
        
        dudxy =  (-(-U(i+2,j+2) + 8*U(i+1,j+2) - 8*U(i-1,j+2) + U(i-2,j+2)) ...
            +8*(-U(i+2,j+1) + 8*U(i+1,j+1) - 8*U(i-1,j+1) + U(i-2,j+1)) ...
            -8*(-U(i+2,j-1) + 8*U(i+1,j-1) - 8*U(i-1,j-1) + U(i-2,j-1)) ...
            + (-U(i+2,j-2) + 8*U(i+1,j-2) - 8*U(i-1,j-2) + U(i-2,j-2)))./(12*hx*12*hy);
        
    end


    function dudx3 = approximateDudx3(fS,i,j,P,hx,hy)
        
        % dudx3 = -dvdx2y
        %       =  dvdy3 - 1/mu*(dvdty + dudy*dvdx + u*dvdxy
        %                       +dvdy*dvdy + v*dvdy2 + dpdxy);
        
        mu    = fS.mu;
        dudy  = fS.dudy;
        dvdy3 = fS.dvdy3;
        dvdy2 = fS.dvdy2;
        dvdx  = fS.dvdx;
        dvdy  = fS.dvdy;
        dvdxy = fS.dvdxy;
        dvdty = fS.dvdty;
        dfydy = fS.dfydy;
        dvdx2y = fS.dvdx2y;
        
        
        dpdy2 = (P(i,j+1) - 2*P(i,j) + P(i,j-1))/(hy^2);
        
        %dudx3 = -dvdx2y(x(i,j),y(i,j),t2);
        dudx3 = dvdy3(x(i,j),y(i,j),t2) - 1/mu*(dvdty(x(i,j),y(i,j),t2) ...
            + dudy(x(i,j),y(i,j),t2).*dvdx(x(i,j),y(i,j),t2) ...
            + u(x(i,j),y(i,j),t2).*dvdxy(x(i,j),y(i,j),t2) ...
            + dvdy(x(i,j),y(i,j),t2).*dvdy(x(i,j),y(i,j),t2) ...
            + v(x(i,j),y(i,j),t2).*dvdy2(x(i,j),y(i,j),t2) ...
            + dpdy2 - dfydy(x(i,j),y(i,j),t2));
%         
       
        
    end


    function U = solveGhostU(fS,BC,Nxg,Nyg,t2,u,v,x,y,hx,hy,dt,dte,dtn,...
            U,...
            VN,VC,VP1,VP2,...
            PN,PC,PP1,PP2,...
            count,tOrder)
        
      dudx2 = fS.dudx2;
      dudx  = fS.dudx;
      dvdx2 = fS.dvdx2;
      
      dvdy2 = fS.dvdy2;
      dvdy  = fS.dvdy;
      dudy2 = fS.dudy2;
      dpdx  = fS.dpdx;
      dpdy  = fS.dpdy;
      
        dudxy = fS.dudxy;
        dvdxy = fS.dvdxy;
        

        
        % boundarys for u
        if BC == 1
            i = 3;
            J = 1:Nyg;
            U(i,J) =u(x(i,J),y(i,J),t2);
       
        
        i = Nxg-2;
        J = 1:Nyg;
        U(i,J) =u(x(i,J),y(i,J),t2);
        end
        
        j  = 3;
        I  = 1:Nxg;
        U(I,j) =u(x(I,j),y(I,j),t2);
        
        j  = Nyg-2;
        I  = 1:Nxg;
        U(I,j) =u(x(I,j),y(I,j),t2);
        
        if BC == 1
        % ghosts in x directions
        i = 1;
        J = 4:Nyg-3;
        iB = i+2;
        
        dir=0;
        
        if count==1
            
            if tOrder==2
                VAppprox = VC;
            elseif tOrder ==4
                VAppprox = 3*VC - 3*VP1 + VP2;
                %VAppprox = VC;
            end
            
            dvdxyApprox = approximateDudxy(iB,J,VAppprox,hx,hy);
            %dvdxyApprox = dvdxy(x(iB,J),y(iB,J),t2);

            U(i,J) = -dvdxyApprox;
            
            if tOrder==2
                PAppprox = PC;
            elseif tOrder ==4
                PAppprox = 3*PC - 3*PP1 + PP2;
            end
            
            dudx3Approx = approximateDudx3(fS,iB,J,PAppprox,hx,hy);
            U(i,J) = dudx3Approx;
        else
            
            dvdxyApprox = approximateDudxy(iB,J,VN,hx,hy);
            %dvdxyApprox = dvdxy(x(iB,J),y(iB,J),t2);
                        
            U(i,J) = -dvdxyApprox;
            
            
            dudx3Approx = approximateDudx3(fS,iB,J,PN,hx,hy);
            U(i,J) = dudx3Approx;
            
        end
        

        
        
        i = Nxg;
        iB = i-2;
        if count==1
            
            if tOrder==2
                VAppprox = VC;
            elseif tOrder ==4
                VAppprox = 3*VC - 3*VP1 + VP2;
                %VAppprox = VC;
            end
            
            dvdxyApprox = approximateDudxy(iB,J,VAppprox,hx,hy);
            %dvdxyApprox = dvdxy(x(iB,J),y(iB,J),t2);
            
            U(i,J) = -dvdxyApprox;
            
        else
            dvdxyApprox = approximateDudxy(iB,J,VN,hx,hy);
            %dvdxyApprox = dvdxy(x(iB,J),y(iB,J),t2);
            U(i,J) = -dvdxyApprox;            
        end
        
        i = 2;
        iB = i+1;
        U(i,J) = -dvdy(x(iB,J),y(iB,J),t2);
        
        i = Nxg-1;
        iB = i-1;
        U(i,J) = -dvdy(x(iB,J),y(iB,J),t2);
        
        end
        
        dir=0;
        % ghosts in y direction
        if BC==1
            I = 4:Nxg-3;
        elseif BC == 6
            I = 3:Nxg-3;
        end
        j = 1;
        U(I,j) =0;
        
        j = Nyg;
        U(I,j) =0;
        
        j = 2;
        jB = j+1;
        if tOrder==2
            dpdxApprox = approximateDpdx(I,jB,count,PN,PC,PP1,PP2,hx,dt,dte,dtn,dir);
        elseif tOrder==4
            dpdxApprox = approximateDpdx(I,jB,count,PN,PC,PP1,PP2,hx,dt,dte,dtn,dir);
           %dpdxApprox = dpdx(x(I,jB),y(I,jB),t2);
        end
        
        if tw == 0
            
            U(I,j) =(1/mu)*( ...
            + dpdxApprox ...  % need to touch this
            - fx(x(I,jB),y(I,jB),t2)) ;
        
        elseif tw == 1
        U(I,j) =(1/mu)*(dudt(x(I,jB),y(I,jB),t2) ...
            +  u(x(I,jB),y(I,jB),t2).*dudx(x(I,jB),y(I,jB),t2) ...
            +  v(x(I,jB),y(I,jB),t2).*dudy(x(I,jB),y(I,jB),t2) ...
            + dpdxApprox ...  % need to touch this
            - fx(x(I,jB),y(I,jB),t2)) ...
            - dudx2(x(I,jB),y(I,jB),t2);
        end
        
        j = Nyg-1;
        jB = j-1;
        if tOrder==2
            dpdxApprox = approximateDpdx(I,jB,count,PN,PC,PP1,PP2,hx,dt,dte,dtn,dir);
        elseif tOrder==4
            dpdxApprox = approximateDpdx(I,jB,count,PN,PC,PP1,PP2,hx,dt,dte,dtn,dir);
           %dpdxApprox = dpdx(x(I,jB),y(I,jB),t2);
        end
        
        if tw == 0
            
            U(I,j) =(1/mu)*( ...
            + dpdxApprox ...  % need to touch this
            - fx(x(I,jB),y(I,jB),t2)) ;
        
        elseif tw == 1
            
        U(I,j) =(1/mu)*(dudt(x(I,jB),y(I,jB),t2) ...
            +  u(x(I,jB),y(I,jB),t2).*dudx(x(I,jB),y(I,jB),t2) ...
            +  v(x(I,jB),y(I,jB),t2).*dudy(x(I,jB),y(I,jB),t2) ...
            + dpdxApprox ...  % need to touch this
            - fx(x(I,jB),y(I,jB),t2)) ...
            - dudx2(x(I,jB),y(I,jB),t2);
        end
        
        if BC == 1
        % corner points
        for i = 1:2
            for j = 1:2
                U(i,j) = u(x(i,j),y(i,j),t2);
            end
        end
        
        % fix corner points for U
        iB = 3;
        jB = 3;
        D1 = hx*dudx(x(iB,jB),y(iB,jB),t2) + hy*dudy(x(iB,jB),y(iB,jB),t2);
        D2 = hx^2/2*dudx2(x(iB,jB),y(iB,jB),t2) + hy^2/2*dudy2(x(iB,jB),y(iB,jB),t2) - hx*hy*dvdy2(x(iB,jB),y(iB,jB),t2);
        U(2,2) = 3/2*D1 + 3*D2;
        
        iB = 3;
        jB = 3;
        iShift = 1;
        jShift = 2;
        D1 = (iShift)*hx*dudx(x(iB,jB),y(iB,jB),t2) + (jShift)*hy*dudy(x(iB,jB),y(iB,jB),t2);
        D2 = ((iShift)*hx)^2/2*dudx2(x(iB,jB),y(iB,jB),t2) + ((jShift)*hy)^2/2*dudy2(x(iB,jB),y(iB,jB),t2) - (iShift)*hx*(jShift)*hy*dvdy2(x(iB,jB),y(iB,jB),t2);
        
        U(iB-iShift,jB-jShift) =  3/2*D1 + 3*D2;
        
        iShift = 2;
        jShift = 1;
        D1 = (iShift)*hx*dudx(x(iB,jB),y(iB,jB),t2) + (jShift)*hy*dudy(x(iB,jB),y(iB,jB),t2);
        D2 = ((iShift)*hx)^2/2*dudx2(x(iB,jB),y(iB,jB),t2) + ((jShift)*hy)^2/2*dudy2(x(iB,jB),y(iB,jB),t2) - (iShift)*hx*(jShift)*hy*dvdy2(x(iB,jB),y(iB,jB),t2);
        
        U(iB-iShift,jB-jShift) =  3/2*D1 + 3*D2;
        
                
        for i = Nxg-1:Nxg
            for j = 1:2
                U(i,j) = u(x(i,j),y(i,j),t2);
            end
        end
        
        for i = 1:2
            for j = Nyg-1:Nyg
                U(i,j) = u(x(i,j),y(i,j),t2);
            end
        end
        
        iB = 3;
        jB = Nyg-2;
        D1 = hx*dudx(x(iB,jB),y(iB,jB),t2) - hy*dudy(x(iB,jB),y(iB,jB),t2);
        D2 = hx^2/2*dudx2(x(iB,jB),y(iB,jB),t2) + hy^2/2*dudy2(x(iB,jB),y(iB,jB),t2) - hx*(-hy)*dvdy2(x(iB,jB),y(iB,jB),t2);
        U(2,Nyg-1) = 3/2*D1 + 3*D2;
        
        for i = Nxg-1:Nxg
            for j = Nyg-1:Nyg
                U(i,j) = u(x(i,j),y(i,j),t2);
            end
        end
        
        end
        
    end



    function V = solveGhostV(V,U,UC,UP1,UP2)
        
        dudx2 = fS.dudx2;
        dudx  = fS.dudx;
        dvdx2 = fS.dvdx2;
        
        dvdy2 = fS.dvdy2;
        dvdy  = fS.dvdy;
        dudy2 = fS.dudy2;
        dpdx  = fS.dpdx;
        dpdy  = fS.dpdy;
        
        dudxy = fS.dudxy;
        dvdxy = fS.dvdxy;
        
        %% boundarys for v
        if BC == 1
            
            i = 3;
            J = 1:Nyg;
            V(i,J) =v(x(i,J),y(i,J),t2);
        
        i = Nxg-2;
        J = 1:Nyg;
        V(i,J) =v(x(i,J),y(i,J),t2);
        end
        
        j  = 3;
        I  = 1:Nxg;
        V(I,j) =v(x(I,j),y(I,j),t2);
        
        j  = Nyg-2;
        I  = 1:Nxg;
        V(I,j) =v(x(I,j),y(I,j),t2);
        
        if BC == 1
        % ghosts in x directions
        i = 1;
        J = 4:Nyg-3;
        iB = i+2;
        V(i,J) =0;
        
        i = Nxg;
        iB = i-2;
        V(i,J) =0;
        
        i = 2;
        iB = i+1;
        dir=1;
        if tOrder==2
            dpdyApprox = approximateDpdx(iB,J,count,PN,PC,PP1,PP2,hy,dt,dte,dtn,dir);
        elseif tOrder==4
            dpdyApprox = approximateDpdx(iB,J,count,PN,PC,PP1,PP2,hy,dt,dte,dtn,dir);
            %            dpdyApprox = dpdy(x(iB,J),y(iB,J),t2);
        end
        
        V(i,J) = (1/mu)*(dvdt(x(iB,J),y(iB,J),t2) ...
            +  u(x(iB,J),y(iB,J),t2).*dvdx(x(iB,J),y(iB,J),t2) ...
            +  v(x(iB,J),y(iB,J),t2).*dvdy(x(iB,J),y(iB,J),t2) ...
            + dpdyApprox ...  % need to touch this
            - fy(x(iB,J),y(iB,J),t2)) ...
            - dvdy2(x(iB,J),y(iB,J),t2);
        
        %V(i,J) = dvdx2(x(iB,J),y(iB,J),t2);
        
        i = Nxg-1;
        iB = i-1;
        if tOrder==2
            dpdyApprox = approximateDpdx(iB,J,count,PN,PC,PP1,PP2,hy,dt,dte,dtn,dir);
        elseif tOrder==4
            dpdyApprox = approximateDpdx(iB,J,count,PN,PC,PP1,PP2,hy,dt,dte,dtn,dir);
            %             dpdyApprox = dpdy(x(iB,J),y(iB,J),t2);
        end
        
        V(i,J) = (1/mu)*(dvdt(x(iB,J),y(iB,J),t2) ...
            +  u(x(iB,J),y(iB,J),t2).*dvdx(x(iB,J),y(iB,J),t2) ...
            +  v(x(iB,J),y(iB,J),t2).*dvdy(x(iB,J),y(iB,J),t2) ...
            + dpdyApprox ...  % need to touch this
            - fy(x(iB,J),y(iB,J),t2)) ...
            - dvdy2(x(iB,J),y(iB,J),t2);
        %V(i,J) =dvdx2(x(iB,J),y(iB,J),t2);
        
        end
        
        % ghosts in y direction
        if BC == 1
            I = 4:Nxg-3;
        elseif BC == 6
            I = 3:Nxg-3;
        end
        j = 1;
        jB = j+2;
        
        if count == 1
            
            if BC == 6 && fS.uvSwitch == 2
                if tOrder==2
                    UAppprox = UC;
                elseif tOrder ==4
                    UAppprox = 3*UC - 3*UP1 + UP2;
                end
                
            else
                UAppprox = U;
            end
            
            if tOrder==2
                dudxyApprox = approximateDudxy(I,jB,UAppprox,hx,hy);
            elseif tOrder==4
                dudxyApprox = approximateDudxy(I,jB,UAppprox,hx,hy);
                %            dudxyApprox = dudxy(x(I,jB),y(I,jB),t2);
            end
        
        else
            dudxyApprox = approximateDudxy(I,jB,U,hx,hy);
        end
        
        V(I,j) =-dudxyApprox; %-dudxy(x(I,jB),y(I,jB),t2);    %% known after calculate u
        
        j = Nyg;
        jB = j-2;
        
        if count == 1
            
            if BC == 6 && fS.uvSwitch == 2
                if tOrder==2
                    UAppprox = UC;
                elseif tOrder ==4
                    UAppprox = 3*UC - 3*UP1 + UP2;
                end
            end
            
            if tOrder==2
                dudxyApprox = approximateDudxy(I,jB,UAppprox,hx,hy);
            elseif tOrder==4
                dudxyApprox = approximateDudxy(I,jB,UAppprox,hx,hy);
                %            dudxyApprox = dudxy(x(I,jB),y(I,jB),t2);
            end
            
        else
            dudxyApprox = approximateDudxy(I,jB,U,hx,hy);
        end
        
        
%         if tOrder==2
%             dudxyApprox = approximateDudxy(I,jB,U,hx,hy);
%         elseif tOrder==4
%             dudxyApprox = approximateDudxy(I,jB,U,hx,hy);
%             %             dudxyApprox = dudxy(x(I,jB),y(I,jB),t2);
%         end
        V(I,j) =-dudxyApprox;%-dudxy(x(I,jB),y(I,jB),t2);
        
        j = 2;
        jB = j+1;
        if tw== 1
        V(I,j) =-dudx(x(I,jB),y(I,jB),t2);
        elseif tw==0 
        V(I,j) = 0;
        end
        
        j = Nyg-1;
        jB = j-1;
        if tw== 1
        V(I,j) =-dudx(x(I,jB),y(I,jB),t2);
        elseif tw==0
            V(I,j) = 0;
        end
        
        if BC == 1
        % corner points
        for i = 1:2
            for j = 1:2
                V(i,j) = v(x(i,j),y(i,j),t2);
            end
        end
        
        
        %% fix corner points for V
        iB = 3;
        jB = 3;
        D1 = hx*dvdx(x(iB,jB),y(iB,jB),t2) + hy*dvdy(x(iB,jB),y(iB,jB),t2);
        D2 = hx^2/2*dvdx2(x(iB,jB),y(iB,jB),t2) + hy^2/2*dvdy2(x(iB,jB),y(iB,jB),t2) - hx*hy*dudx2(x(iB,jB),y(iB,jB),t2);
        V(2,2) = 3/2*D1 + 3*D2;
        
        iB = 3;
        jB = 3;
        iShift = 1;
        jShift = 2;
        D1 = (iShift)*hx*dvdx(x(iB,jB),y(iB,jB),t2) + (jShift)*hy*dvdy(x(iB,jB),y(iB,jB),t2);
        D2 = ((iShift)*hx)^2/2*dvdx2(x(iB,jB),y(iB,jB),t2) + ((jShift)*hy)^2/2*dvdy2(x(iB,jB),y(iB,jB),t2) - (iShift)*hx*(jShift)*hy*dudx2(x(iB,jB),y(iB,jB),t2);
        
        V(iB-iShift,jB-jShift) =  3/2*D1 + 3*D2;
        
        iShift = 2;
        jShift = 1;
        D1 = (iShift)*hx*dvdx(x(iB,jB),y(iB,jB),t2) + (jShift)*hy*dvdy(x(iB,jB),y(iB,jB),t2);
        D2 = ((iShift)*hx)^2/2*dvdx2(x(iB,jB),y(iB,jB),t2) + ((jShift)*hy)^2/2*dvdy2(x(iB,jB),y(iB,jB),t2) - (iShift)*hx*(jShift)*hy*dudx2(x(iB,jB),y(iB,jB),t2);
        
        V(iB-iShift,jB-jShift) =  3/2*D1 + 3*D2;
        
        
        for i = Nxg-1:Nxg
            for j = 1:2
                V(i,j) = v(x(i,j),y(i,j),t2);
            end
        end
        
        for i = 1:2
            for j = Nyg-1:Nyg
                V(i,j) = v(x(i,j),y(i,j),t2);
            end
        end
        
        for i = Nxg-1:Nxg
            for j = Nyg-1:Nyg
                V(i,j) = v(x(i,j),y(i,j),t2);
            end
        end
        end
        
    end

end

