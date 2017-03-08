function P = PressureSolver(t2,fS)

% (t2,beta,g,tref,p,x,y,Nxg,Nyg,...
%     alpha,mubx,muby,hx,hy,fx,fy,...
%     UN,VN,TemN,dubcdt,dvbcdt,dfxdx,dfydy,...
%     Lpl,Lpu,Pp,BC,dtn)

dtn   = fS.dtn;
beta  = fS.beta;
alpha = fS.alpha;

tref = fS.tref;
g    = fS.g;
mubx = fS.mubx;
muby = fS.muby;

x    = fS.x;
y    = fS.y;
Nxg    = fS.Nxg;
Nyg    = fS.Nyg;
hx   = fS.hx;
hy   = fS.hy;

BC   = fS.BC;
if BC == 4
    BC =3;
end

tw = fS.tw;

UN   = fS.UN;
VN   = fS.VN;
TemN = fS.TemN;


p  = fS.p;
fx = fS.fx;
fy = fS.fy;
dubcdt = fS.dubcdt;
dvbcdt = fS.dvbcdt;
dfxdx  = fS.dfxdx;
dfydy  = fS.dfydy;

Lpl = fS.Lpl;
Lpu = fS.Lpu;
Pp  = fS.Pp;


P = zeros(Nxg,Nyg);
RHSp = zeros(Nxg,Nyg);

if BC == 1 || BC == 6
    
    
    dudx2 = fS.dudx2;
    dudx  = fS.dudx;
    dvdy  = fS.dvdy;
    mu    = fS.mu;
    
    if BC == 1
        
    %x-ghost1
    i = 3;
    j = 3:Nyg - 2;
    
    dudxApprox = -dvdy(x(i,j),y(i,j),t2);
    %dudxApprox = (-UN(i+2,j) + 8*UN(i+1,j) - 8*UN(i-1,j) + UN(i-2,j))/(12*hx);
    
    
    dudx2Approx =  -mu*(-(-VN(i+2,j+2) + 8*VN(i+1,j+2) - 8*VN(i-1,j+2) + VN(i-2,j+2)) ...
            +8*(-VN(i+2,j+1) + 8*VN(i+1,j+1) - 8*VN(i-1,j+1) + VN(i-2,j+1)) ...
            -8*(-VN(i+2,j-1) + 8*VN(i+1,j-1) - 8*VN(i-1,j-1) + VN(i-2,j-1)) ...
            + (-VN(i+2,j-2) + 8*VN(i+1,j-2) - 8*VN(i-1,j-2) + VN(i-2,j-2)))./(12*hx*12*hy);
        
    %dudx2Approx = mubx*(-UN(i+2,j) + 16*UN(i+1,j) -30*UN(i,j) + 16*UN(i-1,j)-UN(i-2,j));
    
    
    RHSp(i-1,j) = -dubcdt(x(i,1),y(i,j),t2) ...
        - UN(i,j).*dudxApprox ...
        - VN(i,j).*(-UN(i,j+2) + 8*UN(i,j+1) - 8*UN(i,j-1) + UN(i,j-2))/(12*hy)...
        + dudx2Approx ...
        + muby*(-UN(i,j+2) + 16*UN(i,j+1) -30*UN(i,j) + 16*UN(i,j-1)-UN(i,j-2))...
         ...+ ad21*(-U2(i+2,j) + 16*U2(i+1,j) -30*U2(i,j) + 16*U2(i-1,j)-U2(i-2,j))...
         ...+ ad21*(-U2(i,j+2) + 16*U2(i,j+1) -30*U2(i,j) + 16*U2(i,j-1)-U2(i,j-2))...
        + fx(x(i,1),y(i,j),t2);
    %
    i = Nxg-2;
    j = 3:Nyg - 2;
    
    dudxApprox = -dvdy(x(i,j),y(i,j),t2);
    %dudxApprox = (-UN(i+2,j) + 8*UN(i+1,j) - 8*UN(i-1,j) + UN(i-2,j))/(12*hx);
    
    
    dudx2Approx =  -mu*(-(-VN(i+2,j+2) + 8*VN(i+1,j+2) - 8*VN(i-1,j+2) + VN(i-2,j+2)) ...
        +8*(-VN(i+2,j+1) + 8*VN(i+1,j+1) - 8*VN(i-1,j+1) + VN(i-2,j+1)) ...
        -8*(-VN(i+2,j-1) + 8*VN(i+1,j-1) - 8*VN(i-1,j-1) + VN(i-2,j-1)) ...
        + (-VN(i+2,j-2) + 8*VN(i+1,j-2) - 8*VN(i-1,j-2) + VN(i-2,j-2)))./(12*hx*12*hy);
        
    %dudx2Approx = mubx*(-UN(i+2,j) + 16*UN(i+1,j) -30*UN(i,j) + 16*UN(i-1,j)-UN(i-2,j));
    
    
    
    RHSp(i+1,j) = -dubcdt(x(i,1),y(i,j),t2)...
        -UN(i,j).*dudxApprox...
        - VN(i,j).*(-UN(i,j+2) + 8*UN(i,j+1) - 8*UN(i,j-1) + UN(i,j-2))/(12*hy)...
        + dudx2Approx ...
        + muby*(-UN(i,j+2) + 16*UN(i,j+1) -30*UN(i,j) + 16*UN(i,j-1)-UN(i,j-2))...
        ... + ad21*(-U2(i+2,j) + 16*U2(i+1,j) -30*U2(i,j) + 16*U2(i-1,j)-U2(i-2,j))...
        ... + ad21*(-U2(i,j+2) + 16*U2(i,j+1) -30*U2(i,j) + 16*U2(i,j-1)-U2(i,j-2))...
        + fx(x(i,1),y(i,j),t2);
    end
    
    %y-ghost1
    if BC == 1
        j = 3;
        i = 3:Nxg - 2;
    elseif BC == 6
        j = 3;
        i = 3:Nxg - 3;
    end
    
    if tw == 1
    RHSp(i,j-1) = -dvbcdt(x(i,j),y(1,j),t2)- UN(i,j).*(-VN(i+2,j) + 8*VN(i+1,j) - 8*VN(i-1,j) + VN(i-2,j))/(12*hx)...
        - VN(i,j).*(-VN(i,j+2) + 8*VN(i,j+1) - 8*VN(i,j-1) + VN(i,j-2))/(12*hy) ...
        + mubx*(-VN(i+2,j) + 16*VN(i+1,j) -30*VN(i,j) + 16*VN(i-1,j)- VN(i-2,j)) ...
        + muby*(-VN(i,j+2) + 16*VN(i,j+1) -30*VN(i,j) + 16*VN(i,j-1)- VN(i,j-2))...
        ... + ad21*(-V2(i+2,j) + 16*V2(i+1,j) -30*V2(i,j) + 16*V2(i-1,j)- V2(i-2,j)) ...
        ... + ad21*(-V2(i,j+2) + 16*V2(i,j+1) -30*V2(i,j) + 16*V2(i,j-1)- V2(i,j-2))...
        + fy(x(i,j),y(1,j),t2) ...
        + beta*g*(TemN(i,j)-tref);
    
    elseif tw==0
    RHSp(i,j-1) = - UN(i,j).*(-VN(i+2,j) + 8*VN(i+1,j) - 8*VN(i-1,j) + VN(i-2,j))/(12*hx)...
        - VN(i,j).*(-VN(i,j+2) + 8*VN(i,j+1) - 8*VN(i,j-1) + VN(i,j-2))/(12*hy) ...
        + mubx*(-VN(i+2,j) + 16*VN(i+1,j) -30*VN(i,j) + 16*VN(i-1,j)- VN(i-2,j)) ...
        + muby*(-VN(i,j+2) + 16*VN(i,j+1) -30*VN(i,j) + 16*VN(i,j-1)- VN(i,j-2))...
        ... + ad21*(-V2(i+2,j) + 16*V2(i+1,j) -30*V2(i,j) + 16*V2(i-1,j)- V2(i-2,j)) ...
        ... + ad21*(-V2(i,j+2) + 16*V2(i,j+1) -30*V2(i,j) + 16*V2(i,j-1)- V2(i,j-2))...
        + fy(x(i,j),y(1,j),t2) ...
        + beta*g*(TemN(i,j)-tref);        
    end
    %
    if BC == 1
        j = Nyg-2;
        i = 3:Nxg - 2;
    elseif BC == 6
        j = Nyg-2;
        i = 3:Nxg - 3;
    end

    if tw == 1
    RHSp(i,j+1) = -dvbcdt(x(i,j),y(1,j),t2)- UN(i,j).*(-VN(i+2,j) + 8*VN(i+1,j) - 8*VN(i-1,j) + VN(i-2,j))/(12*hx)...
        - VN(i,j).*(-VN(i,j+2) + 8*VN(i,j+1) - 8*VN(i,j-1) + VN(i,j-2))/(12*hy) ...
        + mubx*(-VN(i+2,j) + 16*VN(i+1,j) -30*VN(i,j) + 16*VN(i-1,j)- VN(i-2,j)) ...
        + muby*(-VN(i,j+2) + 16*VN(i,j+1) -30*VN(i,j) + 16*VN(i,j-1)- VN(i,j-2))...
        ... + ad21*(-V2(i+2,j) + 16*V2(i+1,j) -30*V2(i,j) + 16*V2(i-1,j)- V2(i-2,j)) ...
        ... + ad21*(-V2(i,j+2) + 16*V2(i,j+1) -30*V2(i,j) + 16*V2(i,j-1)- V2(i,j-2))...
        + fy(x(i,j),y(1,j),t2) ...
        + beta*g*(TemN(i,j)-tref);
    
    elseif tw == 0
         RHSp(i,j+1) = - UN(i,j).*(-VN(i+2,j) + 8*VN(i+1,j) - 8*VN(i-1,j) + VN(i-2,j))/(12*hx)...
        - VN(i,j).*(-VN(i,j+2) + 8*VN(i,j+1) - 8*VN(i,j-1) + VN(i,j-2))/(12*hy) ...
        + mubx*(-VN(i+2,j) + 16*VN(i+1,j) -30*VN(i,j) + 16*VN(i-1,j)- VN(i-2,j)) ...
        + muby*(-VN(i,j+2) + 16*VN(i,j+1) -30*VN(i,j) + 16*VN(i,j-1)- VN(i,j-2))...
        ... + ad21*(-V2(i+2,j) + 16*V2(i+1,j) -30*V2(i,j) + 16*V2(i-1,j)- V2(i-2,j)) ...
        ... + ad21*(-V2(i,j+2) + 16*V2(i,j+1) -30*V2(i,j) + 16*V2(i,j-1)- V2(i,j-2))...
        + fy(x(i,j),y(1,j),t2) ...
        + beta*g*(TemN(i,j)-tref);   
        
    end
    
else
    if BC == 2
        
    elseif BC == 3
        
        i = 1:2;
        RHSp(i,:) = p(x(i,:),y(i,:),t2);
        i = Nxg-1:Nxg;
        RHSp(i,:) = p(x(i,:),y(i,:),t2);
        j = 1:2;
        RHSp(:,j) = p(x(:,j),y(:,j),t2);
        j = Nyg-1:Nyg;
        RHSp(:,j) = p(x(:,j),y(:,j),t2);
        
    end
    
end


% interior points and boundaries
if BC == 2
    i = 3:Nxg-3;
    j = 3:Nyg-3;
    
elseif BC == 6
    i = 3:Nxg-3;
    j = 3:Nyg-2;
else
    i = 3:Nxg-2;
    j = 3:Nyg-2;
end
cdx = alpha*(12*mubx);
%cdy = alpha*(12*muby);

if cdx == 0
    cdx = alpha/dtn;%1*0.1/hx^2;
end

% if cdy == 0
%     cdy = alpha/dtn;%1*0.1/hy^2;
% end
if tw == 1
    
forcingTerms =  dfxdx(x(i,j),y(i,j),t2) + dfydy(x(i,j),y(i,j),t2);

elseif tw == 0

forcingTerms = 0;

end

divergence = ((-UN(i+2,j) + 8*UN(i+1,j) - 8*UN(i-1,j) + UN(i-2,j))/(12*hx)...
    + (-VN(i,j+2) + 8*VN(i,j+1) - 8*VN(i,j-1) + VN(i,j-2))/(12*hy));

% divergence = (UN(i+1,j) - UN(i-1,j) )/(2*hx)...
%     + (VN(i,j+1) - VN(i,j-1) )/(2*hy);

RHSp(i,j) = RHSp(i,j) - ((-UN(i+2,j) + 8*UN(i+1,j) - 8*UN(i-1,j) + UN(i-2,j))/(12*hx)).^2 ...
    - 2 * (-UN(i,j+2) + 8*UN(i,j+1) - 8*UN(i,j-1) + UN(i,j-2))...
    .*(-VN(i+2,j) + 8*VN(i+1,j) - 8*VN(i-1,j) + VN(i-2,j))/((12*hx)*(12*hy))...
    - ((-VN(i,j+2) + 8*VN(i,j+1) - 8*VN(i,j-1) + VN(i,j-2))/(12*hy)).^2 ...
    + forcingTerms ...
    + beta*g*(-TemN(i,j+2) + 8*TemN(i,j+1) - 8*TemN(i,j-1) + TemN(i,j-2))/(12*hy) ...
    + cdx*divergence;

if BC == 1
    
    % fill in corner 16 points
    for i=1:2
        for j = 1:2
            RHSp(i,j) = 0;
        end
    end
    for i = Nxg-1:Nxg
        for j = 1:2
            RHSp(i,j) = 0;
        end
    end
    for i = 1:2
        for j = Nyg-1:Nyg
            RHSp(i,j) = 0;
        end
    end
    for i = Nxg-1:Nxg
        for j = Nyg-1:Nyg
            RHSp(i,j) = 0;
        end
    end
end

%from RHS vector
if BC ==2
    
    RHSimp = RHSp(3,3:Nyg-3);
    for i = 4:Nxg-3
        RHSimp = [RHSimp RHSp(i,3:Nyg-3)];
    end
    
else
    
    RHSimp = RHSp(1,1:Nyg);
    for i = 2:Nxg
        RHSimp = [RHSimp RHSp(i,1:Nyg)];
    end
    
end

if BC == 1 || BC == 2 || BC == 6
    RHSimp = [RHSimp 0];
end

if fS.directSolve == 0 
    yt = Lpl\(Pp*RHSimp');
    P2n = Lpu\yt;

elseif fS.directSolve == 1
    Lp = fS.Lp;
    P2n = Lp\(RHSimp');
end



if BC == 2
    for i = 1:Nxg-5
        P(i+2,3:Nyg-3) = P2n((i-1)*(Nyg-5)+1:i*(Nyg-5));
    end
    %ghost points
    P(1,:) = P(Nxg-4,:);
    P(2,:) = P(Nxg-3,:);
    P(Nxg-2,:) = P(3,:);
    P(Nxg-1,:) = P(4,:);
    P(Nxg,:) = P(5,:);
    
    P(:,1) = P(:,Nyg-4);
    P(:,2) = P(:,Nyg-3);
    P(:,Nyg-2) = P(:,3);
    P(:,Nyg-1) = P(:,4);
    P(:,Nyg) = P(:,5);
else
    for i = 1:Nxg
        P(i,:) = P2n((i-1)*Nyg+1:i*Nyg);
    end
end

if BC == 1
    Pt = p(x,y,t2);
    
    %fill the corner points
    i = 1:2;
    j = 1:2;
    P(i,j) =p(x(i,j),y(i,j),t2) + mean(mean(P(3:end-2,3:end-2)-Pt(3:end-2,3:end-2)));
    j = Nyg-1:Nyg;
    P(i,j) =p(x(i,j),y(i,j),t2) + mean(mean(P(3:end-2,3:end-2)-Pt(3:end-2,3:end-2)));
    i = Nxg-1:Nxg;
    j = 1:2;
    P(i,j) =p(x(i,j),y(i,j),t2) + mean(mean(P(3:end-2,3:end-2)-Pt(3:end-2,3:end-2)));
    j = Nyg-1:Nyg;
    P(i,j) =p(x(i,j),y(i,j),t2) + mean(mean(P(3:end-2,3:end-2)-Pt(3:end-2,3:end-2)));
end
end