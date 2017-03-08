function [count,Tem,rhstC] = TempSolver(t1,t2,tem,tembc,...
    ftem,dtemdt,dtemdx,dtemdy,...
    U1,UN,...
    V1,V2,TemC,TemP1,...
    rhstP1,rhstP2,rhstP3,...
    count)
% (t1,t2,tOrder,tem,tembc,x,y,Nxg,Nyg,...
%     hx,hy,al,albx,alby,ftem,dtemdt,dtemdx,dtemdy,...
%     U1,U2,...
%     V1,V2,TemC,TemP1,timeint,BC,bcp,...
%     rhstP1,rhstP2,rhstP3,...
%     count)
%% read in fS

tOrder = fS.tOrder;

al   = fS.al;
albx = fS.albx;
alby = fS.alby;

x    = fS.x;
y    = fS.y;
Nxg    = fS.Nxg;
Nyg    = fS.Nyg;
hx   = fS.hx;
hy   = fS.hy;

BC   = fS.BC; 


timeint = fS.preTime;
if tOrder==2 && count>1
   timeint = fS.corrTime;
end

bcp  = fS.bcp;


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



%% temperature solver
Tem = zeros(Nxg,Nyg);
% boundary points
if BC == 2
    i = 3:Nxg-2;
    j = 3:Nyg-2;
else
    i = 4:Nxg-3;
    j = 4:Nyg-3;
end

rhstC = (ftem(x(i,j),y(i,j),t1) ...
    + albx*(-TemP1(i+2,j) + 16*TemP1(i+1,j) -30*TemP1(i,j) + 16*TemP1(i-1,j)-TemP1(i-2,j)) ...
    + alby*(-TemP1(i,j+2) + 16*TemP1(i,j+1) -30*TemP1(i,j) + 16*TemP1(i,j-1)-TemP1(i,j-2))...
    - U1(i,j).*(-TemP1(i+2,j) + 8*TemP1(i+1,j) - 8*TemP1(i-1,j) + TemP1(i-2,j))/(12*hx)...
    - V1(i,j).*(-TemP1(i,j+2) + 8*TemP1(i,j+1) - 8*TemP1(i,j-1) + TemP1(i,j-2))/(12*hy));

if count == 0
    return
end

if tOrder==2
    
    Tem(i,j) = timeint(1)*rhstC...
        + timeint(2)*rhstP1...
        + TemC(i,j);

elseif tOrder==4

    Tem(i,j) = timeint(1)*rhstC...
    + timeint(2)*rhstP1...
    + timeint(3)*rhstP2...
    + timeint(4)*rhstP3...
    + TemC(i,j);
    
end

switch BC
    case 1
        i = 3;
        Tem(i,:) = tembc(x(i,1),y(i,:),t2);
        i = Nxg-2;
        Tem(i,:) = tembc(x(i,1),y(i,:),t2);
        j = 3;
        Tem(:,j) = tembc(x(:,j),y(1,j),t2);
        j = Nyg-2;
        Tem(:,j) = tembc(x(:,j),y(1,j),t2);
        
        %x-ghost
        for k = 1:2
            i = bcp(1,k);
            j = 4:Nyg - 3;
            index = (-1)^(k-1);
            Temx(:,j) = Tem(i+index*(1:4),j);
            
            rhs =  - dtemdt(x(i,1),y(i,j),t2) - V2(i,j).*(Tem(i,j+1)-Tem(i,j-1))/(2*hy)...
                + al*(Tem(i,j+1) - 2*Tem(i,j) + Tem(i,j-1))/(hy^2) + ftem(x(i,1),y(i,j),t2);
%             Tem(i-index*(1),j) = (rhs - index*U2(i,j).*Temx(1,j)/(2*hx) + al*(Temx(1,j) - 2*Tem(i,j))/(hx^2))...
%                 ./(index*(-U2(i,j))/(2*hx) - al/hx^2);            
            Tem(i-index*(1),j) = (rhs - UN(i,j).*dtemdx(x(i,1),y(i,j),t2) + al*(Temx(1,j) - 2*Tem(i,j))/(hx^2))...
                ./(- al/hx^2);
            
            Tem(i-index*(2),j) = -(-6*Tem(i-index*(1),j)+15*Tem(i,j)-20*Temx(1,j) + 15*Temx(2,j) -6*Temx(3,j) + Temx(4,j));
        end
        
        %y-ghost
        for k = 1:2
            
            j = bcp(2,k);
            i = 4:Nxg - 3;
            index = (-1)^(k-1);
            Temy(i,:) = Tem(i,j+index*(1:4));
            
            rhs =  - dtemdt(x(i,1),y(i,j),t2) - UN(i,j).*(Tem(i+1,j)-Tem(i-1,j))/(2*hx)...
                + al*(Tem(i+1,j) - 2*Tem(i,j) + Tem(i-1,j))/(hx^2) + ftem(x(i,1),y(i,j),t2);
            Tem(i,j-index*(1)) = (rhs - index*V2(i,j).*Temy(i,1)/(2*hy) + al*(Temy(i,1) - 2*Tem(i,j))/(hy^2))...
                ./(-index*V2(i,j)/(2*hy) - al/hy^2);
            
            Tem(i,j-index*(1)) = (rhs - V2(i,j).*dtemdy(x(i,1),y(i,j),t2) + al*(Temy(i,1) - 2*Tem(i,j))/(hy^2))...
                ./(- al/hy^2);
            
            Tem(i,j-index*(2)) = -(-6*Tem(i,j-index*(1))+15*Tem(i,j)-20*Temy(i,1) + 15*Temy(i,2) -6*Temy(i,3) + Temy(i,4));
        end
    case 2
        Tem(bcp(1,1)-1,:) = Tem(bcp(1,2)-1,:);
        Tem(bcp(1,1)-2,:) = Tem(bcp(1,2)-2,:);
        Tem(bcp(1,2)+1,:) = Tem(bcp(1,1)+1,:);
        Tem(bcp(1,2)+2,:) = Tem(bcp(1,1)+2,:);
        
        Tem(:,bcp(2,1)-1) = Tem(:,bcp(2,2)-1);
        Tem(:,bcp(2,1)-2) = Tem(:,bcp(2,2)-2);
        Tem(:,bcp(2,2)+1) = Tem(:,bcp(2,1)+1);
        Tem(:,bcp(2,2)+2) = Tem(:,bcp(2,1)+2);
    case 3
        i = 1:3;
        Tem(i,:) = tembc(x(i,:),y(i,:),t2);
        i = Nxg-2:Nxg;
        Tem(i,:) = tembc(x(i,:),y(i,:),t2);
        j = 1:3;
        Tem(:,j) = tembc(x(:,j),y(:,j),t2);
        j = Nyg-2:Nyg;
        Tem(:,j) = tembc(x(:,j),y(:,j),t2);


end
%fill the corner points
i = 1:2;
j = 1:2;
Tem(i,j) =tem(x(i,j),y(i,j),t2);
j = Nyg-1:Nyg;
Tem(i,j) =tem(x(i,j),y(i,j),t2);
i = Nxg-1:Nxg;
j = 1:2;
Tem(i,j) =tem(x(i,j),y(i,j),t2);
j = Nyg-1:Nyg;
Tem(i,j) =tem(x(i,j),y(i,j),t2);

end