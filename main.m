function errf = main(f,fS)

syms x y t 
u(x,y,t)   = f(1);
v(x,y,t)   = f(2);
p(x,y,t)   = f(3);
tem(x,y,t) = f(4);

CFL       = fS.CFL;
tend      = fS.tEnd;
correct   = fS.numberOfCorrector;
nadapt    = fS.numberOfAdapt;
tOrder    = fS.tOrder;
tExplicit = fS. tExplicit;
tMethod   = fS.tMethod;
extOrder  = fS.extOrder;

if fS.makeMovie == 1
    movieU = VideoWriter('movieU.avi');
    movieV = VideoWriter('movieV.avi');
    movieP = VideoWriter('movieP.avi');
    open(movieU);
    open(movieV);
    open(movieP);
    
            gcf=figure(1);
            %set(gcf,'visible','off');
    
end
%--------------------------------------------------------------------------
%Coefficient specification
%--------------------------------------------------------------------------

al = fS.al;
K  = fS.K;
mu = fS.mu;
tref = fS.tref;
beta = fS.beta;
g = fS.g;

%--------------------------------------------------------------------------
%Twilight zone
%--------------------------------------------------------------------------

%Expression of analytical soln twilight zone method
syms t x y
t = t;
x = x;
y = y;

%--------------------------------------------------------------------------
% time and space discretization
%--------------------------------------------------------------------------
setup(fS);

SymstoFunction

% time
getDt(fS);

%--------------------------------------------------------------------------
%Initial Condition: Interior Region and Boundary
%--------------------------------------------------------------------------
 
dt=fS.dt;
x= fS.x;
y= fS.y;


tP1 = dt*(tOrder-1);
tC  = dt*tOrder;

fS.UP1   = fS.u(x,y,tP1);
fS.VP1   = fS.v(x,y,tP1);
fS.PP1   = fS.p(x,y,tP1);
fS.TemP1 = fS.tem(x,y,tP1);

fS.UC   = fS.u(x,y,tC);
fS.VC   = fS.v(x,y,tC);
fS.PC   = fS.p(x,y,tC);
fS.TemC = fS.tem(x,y,tC);
        
textr = tP1-dt;

tP2 = dt*(tOrder-2);
tP3 = dt*(tOrder-3);

fS.UP2   = fS.u(x,y,tP2);
fS.VP2   = fS.v(x,y,tP2);
fS.PP2   = fS.p(x,y,tP2);
fS.TemP2 = fS.tem(x,y,tP2);

fS.UP3   = fS.u(x,y,tP3);
fS.VP3   = fS.v(x,y,tP3);
fS.PP3   = fS.p(x,y,tP3);
fS.TemP3 = fS.tem(x,y,tP3);

if tOrder==4
    
   textr = tP3-dt;
end    
    
fS.Pextr = fS.p(x,y,textr);
    

if tOrder==2
    
    fS.PP2 = fS.Pextr; 
    fS.PP3 = 0;
end
  

BC = fS.BC;

Nxg = fS.Nxg;
Nyg = fS.Nyg;
hx  = fS.hx;
hy  = fS.hy;
    L =  formLHS(fS);
if (BC == 1)  || (BC == 6)
    Lp = formLHSp1;
elseif BC == 2
    Lp = formLHSp2;
elseif (BC == 3) || (BC == 4)
    Lp = formLHSp3;    
end

fS.Lp = Lp;

if fS.directSolve == 0
    tic
    [Lpl,Lpu,Pp] = lu(Lp);
    fS.Lpl = Lpl;
    fS.Lpu = Lpu;
    fS.Pp  = Pp;
    toc
end

if  tExplicit == 0
    
    imTime    = fS.imTime;
    
    Lp = formLHSVelU1(fS);
    [Lul,Luu,Up] = lu(Lp);
    
    fS.Lu = Lp;
    
    Lp = formLHSVelV1(fS);
    [Lvl,Lvu,Vp] = lu(Lp);
    
    fS.Lv = Lp;

    if fS.directSolve == 0
        
        fS.Lul = Lul;
        fS.Luu = Luu;
        fS.Up  = Up;
        
        fS.Lvl = Lvl;
        fS.Lvu = Lvu;
        fS.Vp  = Vp;
        
    end
   
end
    

%Initialization
fS.rhsuP1= zeros(Nxg,Nyg);
fS.rhsvP1= zeros(Nxg,Nyg);
fS.rhstP1= zeros(Nxg,Nyg);

%if tOrder==4
fS.rhsuP2= zeros(Nxg,Nyg);
fS.rhsvP2= zeros(Nxg,Nyg);
fS.rhstP2= zeros(Nxg,Nyg);
fS.rhsuP3= zeros(Nxg,Nyg);
fS.rhsvP3= zeros(Nxg,Nyg);
fS.rhstP3= zeros(Nxg,Nyg);

%end


%% initilize the rhs for previous timesteps
count = 0;

fS.UN = fS.UP1;
fS.VN = fS.VP1;
fS.PN = fS.PP1;
fS.TemN = fS.TemP1;

[~,~,~,rhsuP1,rhsvP1] = VelocitySolver(tP1,tP1,count,fS);

fS.rhsuP1=rhsuP1;
fS.rhsvP1=rhsvP1;

if tOrder==4
    
    fS.UN = fS.UP2;
    fS.VN = fS.VP2;
    fS.PN = fS.PP2;
    fS.TemN = fS.TemP2;
    [~,~,~,rhsuP2,rhsvP2] = VelocitySolver(tP2,tP2,count,fS);
    fS.rhsuP2=rhsuP2;
    fS.rhsvP2=rhsvP2;

    fS.UN = fS.UP3;
    fS.VN = fS.VP3;
    fS.PN = fS.PP3;
    fS.TemN = fS.TemP3;
    [~,~,~,rhsuP3,rhsvP3] = VelocitySolver(tP3,tP3,count,fS);
    fS.rhsuP3=rhsuP3;
    fS.rhsvP3=rhsvP3;
    
end
%--------------------------------------------------------------------------
% Start time marching
%--------------------------------------------------------------------------
tic
t1 = tC-dt;
t = [];
dtn = dt; 
dte = dt;
n = 1;
t2= 0;

% if fS.makeMovie == 1
%     surf(UN);
%     axis tight manual
%     set(gca,'nextplot','replacechildren'); 
% end


while (t2 - tend) < eps-dt/4
    n = n+1;

    count = 1;
    t1 = t1 + dt;
    t2 = t1 + dtn;%the time of stage2 in each time step

    t = [t t2];

    %% Predictor
    fS.UN = fS.UC;
    fS.VN = fS.VC;
    fS.PN = fS.PC;
    
    [count,U2,V2,rhsuC,rhsvC,maxgrad] = VelocitySolver(t1,t2,count,fS);
    
    %     [count,Tem2,rhst0] = TempSolver(t1,t2,tem,tembc,x,y,Nxg,Nyg,...
    %         hx,hy,al,albx,alby,ftem,dtemdt,dtemdx,dtemdy,...
    %         U1,U2,...
    %         V1,V2,Tem1,Tem1,ab2,BC,bcp,rhst0,count);
    %
    
    %U2 = fS.u(x,y,t2);
    %V2 = fS.v(x,y,t2);
    Tem2 = fS.tem(x,y,t2);
    
    fS.UN   = U2;
    fS.VN   = V2;
    fS.TemN = Tem2;
    
    %if correct>0
        if tOrder==2
            fS.rhsuP1 = rhsuC;
            fS.rhsvP1 = rhsvC;
        elseif tOrder ==4
            fS.rhsuP3 = fS.rhsuP2;
            fS.rhsuP2 = fS.rhsuP1;
            fS.rhsuP1 = rhsuC;
            
            fS.rhsvP3 = fS.rhsvP2;
            fS.rhsvP2 = fS.rhsvP1;
            fS.rhsvP1 = rhsvC;
        end
    %end
    
    P2 = PressureSolver(t2,fS);
    %P2 = fS.p(x,y,t2);
    
    fS.PN   = P2;
    
    for k = 1:correct
        
        [count,U2,V2] = VelocitySolver(t2,t2,count,fS);
        
        Tem2 = fS.tem(x,y,t2);
        
        fS.UN   = U2;
        fS.VN   = V2;
        fS.TemN = Tem2;
        
        P2 = PressureSolver(t2,fS);
        %P2 = fS.p(x,y,t2);
        
        fS.PN   = P2;
    end
    %
    dte = dt;
    dt = dtn;
    
    %% update solution
    
    if tOrder==2
        
        fS.PP2 = fS.PP1;
        fS.PP1 = fS.PC;
        
        fS.UC   = U2;

        fS.VP2 = fS.VP1;
        fS.VP1 = fS.VC;
        fS.VC   = V2;
        
        fS.PC   = P2;
        fS.TemC = Tem2;
        
        
    elseif tOrder==4
        
        fS.Pextr = fS.PP3;
        fS.PP3   = fS.PP2;
        fS.PP2   = fS.PP1;
        fS.PP1   = fS.PC;
        
        fS.UP3 = fS.UP2;
        fS.UP2 = fS.UP1;
        fS.UP1 = fS.UC;
        fS.UC  = U2;
        
        fS.VP3 = fS.VP2;
        fS.VP2 = fS.VP1;
        fS.VP1 = fS.VC;
        fS.VC  = V2;
        
        
        
        fS.PC = P2;
        fS.TemC = Tem2;
              
    end
    
    %% update time step $ time integration scheme
    if CFL ~= 0
        if mod(n,nadapt) == 1 || nadapt == 1
            ilambda = max(max(U2))*(1+2/3)/hx + max(max(V2))*(1+2/3)/hy;
            if max(abs(ad41),abs(ad42)) > 0
                rlambda = 4*max(mu,al)*(1+1/3)/hx^2 + 4*max(mu,al)*(1+1/3)/hy^2 + 2*(ad41 + ad42*maxgrad)*16;
            else
                rlambda = 4*max(mu,al)*(1+1/3)/hx^2 + 4*max(mu,al)*(1+1/3)/hy^2;
            end
            dtn = sqrt(CFL/((ilambda/1.28)^2 + (rlambda/1.6)^2));
            tntemp = ceil((tend - t2)/dtn);
            dtn = (tend - t2)/tntemp;
            
            preTime = [dtn + dtn^2/dt/2, -dtn^2/dt/2];
            corrTime = [dtn/2,dtn/2];
            
        end
        if mod(n,nadapt) == 2
            preTime = [dtn + dtn^2/dt/2, -dtn^2/dt/2];
        end
    end
    
    
    
    %% difference between ture soln and approxn
    Ut = fS.u(x,y,t2);
    Vt = fS.v(x,y,t2);
    Pt = fS.p(x,y,t2);
    Tt = fS.tem(x,y,t2);
    
    err(1,n+1) = max(max(abs(U2(3:end-2,3:end-2)-Ut(3:end-2,3:end-2))));
    err(2,n+1) = max(max(abs(V2(3:end-2,3:end-2)-Vt(3:end-2,3:end-2))));
    err(3,n+1) = max(max(abs(Tem2(3:end-2,3:end-2)-Tt(3:end-2,3:end-2))));
    err(4,n+1) = max(max(abs(P2(3:end-2,3:end-2)-Pt(3:end-2,3:end-2))...
        -abs(mean(mean(P2(3:end-2,3:end-2)-Pt(3:end-2,3:end-2))))));
    i = 3:Nxg-3;
    j = 3:Nyg-3;
    
    div = abs((-U2(i+2,j) + 8*U2(i+1,j) - 8*U2(i-1,j) + U2(i-2,j))/(12*hx)...
        + (-V2(i,j+2) + 8*V2(i,j+1) - 8*V2(i,j-1) + V2(i,j-2))/(12*hy));
    diverr(n+1) = max(max(div));
    
    vorticity = abs((-V2(i+2,j) + 8*V2(i+1,j) - 8*V2(i-1,j) + V2(i-2,j))/(12*hx)...
        - (-U2(i,j+2) + 8*U2(i,j+1) - 8*U2(i,j-1) + U2(i,j-2))/(12*hy));
    
    momentumerr(1,n+1) = sum(sum(U2(i,j)))-sum(sum(Ut(i,j)));
    momentumerr(2,n+1) = sum(sum(V2(i,j)))-sum(sum(Vt(i,j)));
    
    energyerr(1,n+1) = sum(sum(U2(i,j).^2))-sum(sum(Ut(i,j).^2));
    energyerr(2,n+1) = sum(sum(V2(i,j).^2))-sum(sum(Vt(i,j).^2));
    
    peakuv(1,n+1)= max(max(U2(i,j)));
    peakuv(2,n+1)= max(max(V2(i,j)));
    
    
    %if mod(n,20) == 0
    errorcalculation
    %end
    
    toc
if fS.plotting == 1    
    if abs(t2 - .5*tend) <= dtn
        pl = figure(1);
        debugPlot
        print(pl, '-depsc',[path 't12err' '.eps']);
        
        pl = figure(2);
        contourf(div)
        colorbar
        print(pl, '-depsc',[path 't12div' '.eps']);        
    end
    
    if abs(t2 - .75*tend) <= dtn
        pl = figure(1);
        debugPlot
        print(pl, '-depsc',[path 't34err' '.eps']);
        
        pl = figure(2);
        contourf(div)
        colorbar
        print(pl, '-depsc',[path 't34div' '.eps']);
    end
end    
    
%    if abs(t2 - 4.99) <= dtn
%        pl = figure(1);
%        debugPlot
        
%        pl = figure(2);
%        contourf(div)
%        colorbar
%    end
if fS.plotting == 2
 %debugPlot
 figure(1)
 %surf(x(3:end-2,3:end-2),y(3:end-2,3:end-2),U2(3:end-2,3:end-2)-Ut(3:end-2,3:end-2));
 surf(x,y,U2-Ut);
 figure(2)
 % surf(x(3:end-2,3:end-2),y(3:end-2,3:end-2),V2(3:end-2,3:end-2)-Vt(3:end-2,3:end-2));
  surf(x,y,V2-Vt);
 figure(3)
  surf(x,y,P2-Pt);
 pause
end

if fS.makeMovie == 1
    
    if mod(n,1)==0;
        
        %surf(x(3:end-2,3:end-2),y(3:end-2,3:end-2),U2(3:end-2,3:end-2)); 
        surf(x(3:Nxg-3,3:Nyg-3),y(3:Nxg-3,3:Nyg-3),vorticity);    
        view(0,90);
        shading interp;
        frame = getframe(gcf);
        writeVideo(movieU,frame);
        

        surf(x(3:end-2,3:end-2),y(3:end-2,3:end-2),V2(3:end-2,3:end-2))
        view(0,90);
        shading interp;
        frame = getframe(gcf);
        writeVideo(movieV,frame);
        
        surf(x(3:end-2,3:end-2),y(3:end-2,3:end-2),P2(3:end-2,3:end-2))
        view(0,90);
        shading interp;
        frame = getframe(gcf);
        writeVideo(movieP,frame);
        
    end
end

    if max(err(1:4,n+1)) > 1e+5
        break
    end
    
%end
end

    function SymstoFunction
        
        if fS.tw == 0
            
            if fS.twilightZone == 7
                G  = 1;
            elseif fS.twilightZone == 8
                G  = 0;
            end
            
            fx = G + eps*x + eps*y + eps*t;
            fy = eps*x + eps*y + + eps*t;
            
            if fS.twilightZone == 8
                G = 1;
                H = fS.domain(2,2);
                
                fS.u = @(x,y,t) eps*x + eps*y + eps*t+ ((y>0).*G*H^2/(2).*(1 - y.^2/H^2)) ...
                    + ((y<=0).*(-G)*H^2/(2).*(1 - y.^2/H^2));

                if fS.BC == 2
                    
                    fS.u = @(x,y,t) eps*x + eps*y + eps*t+ ((y>0).*G*H^2/(2).*(1 )) ...
                        + ((y<=0).*(-G)*H^2/(2).*(1 ));
                    
                end
                
            else
                fS.u  = matlabFunction(u);
            end
            
            fS.v  = matlabFunction(v);
            fS.p  = matlabFunction(p);
            fS.tem = matlabFunction(tem);
            
            fS.fx  = matlabFunction(fx);
            fS.fy  = matlabFunction(fy);
            
            fS.fxE  = matlabFunction(fx);
            fS.fyE  = matlabFunction(fy);            
            
        return 
        end
        
        dtemdt = diff(tem,t);
        dtemdx = diff(tem,x);
        dtemdy = diff(tem,y);
        
        dtemdx2= diff(tem,x,2);
        dtemdy2= diff(tem,y,2);
        ftem = dtemdt + u*dtemdx + v*dtemdy - al*dtemdx2 - al*dtemdy2;
        
        dudt = diff(u,t);
        dvdt = diff(v,t);
        
        
        dudx = diff(u,x);
        dudxt = diff(dudx,t);
        
        dudx2= diff(u,x,2);
        dudx3 = diff(u,x,3);
        dudxy = diff(dudx,y);
        dudy3 = diff(u,y,3);
        
        dvdx = diff(v,x);
        dvdx2= diff(v,x,2);
        
        dudy = diff(u,y);
        dudy2= diff(u,y,2);
        
        dvdy = diff(v,y);
        dvdy2= diff(v,y,2);
        dvdy3= diff(v,y,3);
        dvdty = diff(dvdy,t);
        dvdx2y = diff(dvdx2,y);
        
        
        
        dpdx = diff(p,x);
        dpdx2= diff(p,x,2);
        dpdy = diff(p,y);
        dpdy2 = diff(p,y,2);
        
        
        fx = dudt + u*dudx +v*dudy + dpdx - mu*dudx2 - mu*dudy2;
        fy = dvdt + u*dvdx +v*dvdy + dpdy - mu*dvdx2 - mu*dvdy2 - beta*g*(tem-tref);
        %pause
        dfxdx = diff(fx,x);
        dfydy = diff(fy,y);
        
        fxE = dudt + u*dudx +v*dudy + dpdx;
        fyE = dvdt + u*dvdx +v*dvdy + dpdy - beta*g*(tem-tref);
        
        fxI = - mu*dudx2 - mu*dudy2;
        fyI = - mu*dvdx2 - mu*dvdy2;
            
        % Expression for the boundary condition
        ubc = u;
        vbc = v;
        dubcdt = diff(ubc,t);
        dvbcdt = diff(vbc,t);
        
        dvdxy = diff(dvdx,y);
        
        %%mu dudx3
        %mu*dudx3 + dvdyt - dvdy*dvdy + u*dvdxy - dvdx*dudy + v*dvdy2 - dpdx2 - mu*dvdy3 + dfxdx
        
        %-dudxt + dudy*dvdx - u*dudx2  + dudx*dudx - v*dudxy + dpdy2 + mu*dudx3 - mu*dvdy3 - dfydy - g*beta*dtemdy
        
        % functionalized the expressions
      
        fS.u  = matlabFunction(u);
        fS.v  = matlabFunction(v);
        fS.p  = matlabFunction(p);
        fS.tem = matlabFunction(tem);
        
%         ubc = u;
%         vbc = v;
%         tembc = tem;
        
        fS.fx    = matlabFunction(fx);
        fS.fy    = matlabFunction(fy);
        fS.fxE   = matlabFunction(fxE);
        fS.fyE   = matlabFunction(fyE);
        fS.fxI   = matlabFunction(fxI);
        fS.fyI   = matlabFunction(fyI);
        
%        ftem    = matlabFunction(ftem);
        
        fS.dfxdx    = matlabFunction(dfxdx);
        fS.dfydy    = matlabFunction(dfydy);
        fS.dubcdt   = matlabFunction(dubcdt);
        fS.dvbcdt   = matlabFunction(dvbcdt);
        
        fS.dudt = matlabFunction(dudt);
        fS.dvdt = matlabFunction(dvdt);
        
        %dvdyt = matlabFunction(dvdyt);
        
        fS.dudx  = matlabFunction(dudx);
        fS.dudy  = matlabFunction(dudy);
        fS.dvdx  = matlabFunction(dvdx);
        fS.dvdy  = matlabFunction(dvdy);
        fS.dudxy = matlabFunction(dudxy);
        fS.dvdxy = matlabFunction(dvdxy);
        fS.dvdty = matlabFunction(dvdty);
        fS.dvdx2y = matlabFunction(dvdx2y);
        
        fS.dudx2 = matlabFunction(dudx2);
        fS.dudy2 = matlabFunction(dudy2);
        fS.dvdx2 = matlabFunction(dvdx2);
        fS.dvdy2 = matlabFunction(dvdy2);
        fS.dvdy3 = matlabFunction(dvdy3);

        fS.dpdx  = matlabFunction(dpdx);

        fS.dpdy  = matlabFunction(dpdy);
        
%         
%         dudxt = matlabFunction(dudxt);
%         dudx2 = matlabFunction(dudx2);
%         dvdy3 = matlabFunction(dvdy3);
%         dpdx2 = matlabFunction(dpdx2);
%         dudxy = matlabFunction(dudxy);
%         
%         dudx3 = matlabFunction(dudx3);
%         
%         dtemdt = matlabFunction(dtemdt);
%         dtemdx = matlabFunction(dtemdx);
%         dtemdy = matlabFunction(dtemdy);

        
        
    end
    function [L1] = formLHSp1

        M  = Nyg*Nxg;
        
        L1 = spalloc(M,M,7*M);
        
        if BC == 1
            
            for ix = 1:2
                
                L1(ix,ix) = 1;
            end
            
            for ix = 3:Nyg-2
                
                L1(ix,ix        ) =  1;
                L1(ix,ix +   Nyg) = -4;
                L1(ix,ix + 2*Nyg) =  6;
                L1(ix,ix + 3*Nyg) = -4;
                L1(ix,ix + 4*Nyg) =  1;
                
            end
            
            for ix = Nyg-1:Nyg+2
                
                L1(ix,ix) = 1;
                
            end
            
            for ix = Nyg + 3: Nyg + Nyg-2
                
                L1(ix,ix        ) =  (-8)/(12*hx);
                L1(ix,ix -   Nyg) =  ( 1)/(12*hx);
                L1(ix,ix + 2*Nyg) =  ( 8)/(12*hx);
                L1(ix,ix + 3*Nyg) =  (-1)/(12*hx);
                
            end
            
            for ix = 2*Nyg-1:2*Nyg
                
                L1(ix,ix) = 1;
                
            end
            
            
        elseif BC ==6
            
            for ix = 1:2*Nyg;
                
                L1(ix,ix        )     =   1;
                L1(ix,ix + M - 5*Nyg) =  -1;
                
            end
            
        end
        
        
        if BC == 1
            
            indexInterior = M-2*Nyg;
        
        elseif BC == 6
            
            indexInterior = M-3*Nyg;
            
        end
        
        for ix = 2*Nyg+1:indexInterior
            
            if mod(ix - 2*Nyg,Nyg) == 1
                
                L1(ix,ix    ) =  1;
                L1(ix,ix + 1) = -4;
                L1(ix,ix + 2) =  6;
                L1(ix,ix + 3) = -4;
                L1(ix,ix + 4) =  1;                
                
            elseif mod(ix - 2*Nyg,Nyg) == 0
                
                L1(ix,ix    ) =  1;
                L1(ix,ix - 1) = -4;
                L1(ix,ix - 2) =  6;
                L1(ix,ix - 3) = -4;
                L1(ix,ix - 4) =  1;
                
            elseif mod(ix - 2*Nyg,Nyg) ==  2
                
                L1(ix,ix    ) =  (-8)/(12*hy);
                L1(ix,ix - 1) =  ( 1)/(12*hy);
                L1(ix,ix + 2) =  ( 8)/(12*hy);
                L1(ix,ix + 3) =  (-1)/(12*hy);
                
            elseif mod(ix - 2*Nyg,Nyg) ==  Nyg-1
                
                L1(ix,ix    ) =  ( 8)/(12*hy);
                L1(ix,ix + 1) =  (-1)/(12*hy);
                L1(ix,ix - 2) =  (-8)/(12*hy);
                L1(ix,ix - 3) =  ( 1)/(12*hy);
                
                
            else
                L1(ix,ix   ) =  -(30/(12*hy^2)+30/(12*hx^2));
                L1(ix,ix +1) =  16/(12*hy^2);
                L1(ix,ix -1) =  16/(12*hy^2);
                L1(ix,ix +2) =    -1/(12*hy^2);
                L1(ix,ix -2) =    -1/(12*hy^2);
                
                L1(ix,ix + Nyg) =   16/(12*hx^2);
                L1(ix,ix - Nyg) =   16/(12*hx^2);
                L1(ix,ix + 2*Nyg) =   -1/(12*hx^2);
                L1(ix,ix - 2*Nyg) =   -1/(12*hx^2);
                
            end
            
        end
        
        
        if BC == 1
            
            for ix = M-2*Nyg + 1 : M-2*Nyg + 2
                L1(ix,ix) = 1;
            end
            
            for ix = M-2*Nyg + 3 : M-Nyg-2
                L1(ix,ix        ) =  ( 8)/(12*hx);
                L1(ix,ix +   Nyg) =  (-1)/(12*hx);
                L1(ix,ix - 2*Nyg) =  (-8)/(12*hx);
                L1(ix,ix - 3*Nyg) =  ( 1)/(12*hx);
            end
            
            for ix = M-Nyg-1 : M - Nyg + 2
                L1(ix,ix)=  1;
            end
            
            
            for ix = M - Nyg + 3: M-2
                L1(ix,ix)     =  1;
                L1(ix,ix -   Nyg) = -4;
                L1(ix,ix - 2*Nyg) =  6;
                L1(ix,ix - 3*Nyg) = -4;
                L1(ix,ix - 4*Nyg) =  1;
            end
            
            for ix = M-1:M
                L1(ix,ix)=1;
            end
            
        elseif BC ==6
            
            for ix = (M-3*Nyg + 1):M;
                
                L1(ix,ix        )     =   1;
                L1(ix,ix - M + 5*Nyg) =  -1;
                
            end
            
        end
        
        
        %--------------------------------------------------------------------------
        %The compatability condition for pressure equation
        %--------------------------------------------------------------------------
        
        comp = ones(1,M);
        L1 = [L1;comp];
        comp = ones(1,M+1)';
        L1 = [L1 comp];
        L1(end,end) = 0;
        
        L1 = sparse(L1);
        
    end

 function [L1] = formLHSp2
     
        dim = [Nxg-5,Nyg-5];
        M = dim(1)*dim(2);
        Ix = speye(dim(1));
        Iy = speye(dim(2));
        e1 = ones(dim(1),1);
        e2 = ones(dim(2),1);
        
        D1x = spdiags([-1/(12*hx^2)*e1 16/(12*hx^2)*e1 -30/(12*hx^2)*e1 16/(12*hx^2)*e1 -1/(12*hx^2)*e1]...
            , [-2 -1 0 1 2], dim(1),dim(1));
        D1x(1,dim(1)) = 16/(12*hx^2);
        D1x(1,dim(1)-1) = -1/(12*hx^2);
        D1x(2,dim(1)) = -1/(12*hx^2);

        D1x(dim(1),1) = 16/(12*hx^2);
        D1x(dim(1),2) = -1/(12*hx^2);
        D1x(dim(1)-1,1) = -1/(12*hx^2);
        
        D1y = spdiags([-1/(12*hy^2)*e2 16/(12*hy^2)*e2 -30/(12*hy^2)*e2 16/(12*hy^2)*e2 -1/(12*hy^2)*e2]...
            , [-2 -1 0 1 2], dim(2),dim(2));
        
        D1y(1,dim(2)) = 16/(12*hy^2);
        D1y(1,dim(2)-1) = -1/(12*hy^2);
        D1y(2,dim(2)) = -1/(12*hy^2);

        D1y(dim(2),1) = 16/(12*hy^2);
        D1y(dim(2),2) = -1/(12*hy^2);
        D1y(dim(2)-1,1) = -1/(12*hy^2);

        L1 = kron(D1x,Iy) + kron(Ix,D1y);        
        %--------------------------------------------------------------------------
        %The compatability condition for pressure equation
        %--------------------------------------------------------------------------
        
        comp = ones(1,M);
        L1 = [L1;comp];
        comp = ones(1,M+1)';
        L1 = [L1 comp];
        L1(end,end) = 0;
    end

    function [L1] = formLHSp3
        M  = (Nyg)*(Nxg);
        %L1 = zeros(M,M);
                L1 = spalloc(M,M,7*M);
        
        for ix = 1:2*Nyg
            L1(ix,ix) = 1;
        end
        
        for ix = 2*Nyg+1:M-2*Nyg
            if mod(ix - 2*Nyg,Nyg) == 1
                L1(ix,ix   ) = 1;
            elseif mod(ix - 2*Nyg,Nyg) == 2
                L1(ix,ix   ) = 1;
            elseif mod(ix - 2*Nyg,Nyg) == 0
                L1(ix,ix   ) = 1;
            elseif mod(ix - 2*Nyg,Nyg) == Nyg-1
                L1(ix,ix   ) = 1;
            else
                L1(ix,ix   ) =  -(30/(12*hy^2)+30/(12*hx^2));
                
                L1(ix,ix + 1) =  16/(12*hy^2);
                L1(ix,ix - 1) =  16/(12*hy^2);
                L1(ix,ix + 2) =    -1/(12*hy^2);
                L1(ix,ix - 2) =    -1/(12*hy^2);
                
                L1(ix,ix + Nyg) =   16/(12*hx^2);
                L1(ix,ix - Nyg) =   16/(12*hx^2);
                L1(ix,ix + 2*Nyg) =   -1/(12*hx^2);
                L1(ix,ix - 2*Nyg) =   -1/(12*hx^2);
            end
            
        end
        
        for ix = M-2*Nyg+1:M
            L1(ix,ix) = 1;
        end

    end


    function errorcalculation
        fprintf('at time %8.5f \n',t2)
        fprintf('l2 of U = %8.3e \n',err(1,n+1))
        fprintf('l2 of V = %8.3e \n',err(2,n+1))
        fprintf('l2 of T = %8.3e \n',err(3,n+1))
        fprintf('l2 of P = %8.3e \n',err(4,n+1))
        
    end

%if plotting == 1
    close all
    %p = figure('Visible','off');
    pl = figure(1);
    %t = linspace(dt,tend,tn-1);
    label = 'Error vs t';
    plot(t,log10(err(1,3:end)),'LineWidth',2)
    hold on
    plot(t,log10(err(2,3:end)),'r--','LineWidth',2)
    %plot(t,log10(err(3,3:end)),'k','LineWidth',2)
    plot(t,log10(err(4,3:end)),'g','LineWidth',2)
    plot(t,log10(abs(diverr(3:end))),'m','LineWidth',2)

    set(gca, 'FontSize', 20)
    xlabel('t');
    ylabel('log10(Error)');
    title(label,'FontSize',20);
    legend('log10(|U-Uexact|_\infty)','log10(|V-Vexact|_\infty)','log10(|P-Pexact|_\infty)','log10(|\nabla \cdot u|_\infty)','location','southeast')
    print(pl, '-dpng',[fS.path '.png']);
%     end

    function debugPlot        
        Tt = fS.tem(x,y,t2);
        Ut = fS.u(x,y,t2);
        Vt = fS.v(x,y,t2);
        Pt = fS.p(x,y,t2);
        subplot(2,2,1);
        surf(x,y,V2-Vt)
        shading interp
        title('error V','FontSize',16);
        subplot(2,2,2);
        surf(x,y,U2-Ut)
        shading interp
        title('error U','FontSize',16);
%         subplot(2,2,3);
%         surf(x,y,Tem2-Tt)
%         title('error T','FontSize',16);
        subplot(2,2,3);
        surf(x,y,abs(P2-Pt) - abs(mean(mean(P2-Pt))))
        shading interp
        title('error P','FontSize',16);        
    end

pl=figure(2);
debugPlot
print(pl, '-dpng',[fS.path 'Error' '.png']); 

% pl=figure(3);
% plot(t(:),energyerr(2,3:end)+energyerr(1,3:end),'g','LineWidth',2)
% hold on
% plot(t(:),energyerr(1,3:end),'b','LineWidth',2)
% plot(t(:),energyerr(2,3:end),'r','LineWidth',2)
% set(gca, 'FontSize', 20)
% xlabel('t');
% ylabel('Energy - Energy_{exact}');
% title(label,'FontSize',20);
% legend('U^2+V^2','U^2','V^2')
% print(pl, '-depsc',[fS.path 'energy' '.eps']);
% 
% pl=figure(4);
% plot(t(:),momentumerr(1,3:end),'b','LineWidth',2)
% hold on
% plot(t(:),momentumerr(2,3:end),'r','LineWidth',2)
% set(gca, 'FontSize', 20)
% xlabel('t');
% ylabel('momentum - momentum_{exact}');
% title(label,'FontSize',20);
% legend('U','V')
% print(pl, '-depsc',[fS.path 'momentum' '.eps']);
% 
% 
% pl=figure(5);
% plot(t(:),peakuv(1,3:end),'b','LineWidth',2)
% hold on
% plot(t(:),peakuv(2,3:end),'r','LineWidth',2)
% set(gca, 'FontSize', 20)
% xlabel('t');
% ylabel('peak - peak_{exact}');
% title(label,'FontSize',20);
% legend('U','V')
% print(pl, '-depsc',[fS.path 'peak' '.eps']);


close(movieU);
close(movieV);
close(movieP);

%no temperature
for i = 1:2
    errf(i) = err(i,end);
end
errf(3)=err(4,end);

end
