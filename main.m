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
%extOrder  = fS.extOrder;
ad41 = fS.ad41;
ad42 = fS.ad42;

if fS.makeMovie == 1
    movieU = VideoWriter('movieU.avi');
    movieV = VideoWriter('movieV.avi');
    movieP = VideoWriter('movieP.avi');
    open(movieU);
    open(movieV);
    open(movieP);
    
    %gcf=figure(1);
    gcf = figure('Visible','off');
            %set(gcf,'visible','off');
    
end
%--------------------------------------------------------------------------
%Coefficient specification
%--------------------------------------------------------------------------

al = fS.al;
% K  = fS.K;
mu = fS.mu;
tref = fS.tref;
beta = fS.beta;
g = fS.g;

%--------------------------------------------------------------------------
%Twilight zone
%--------------------------------------------------------------------------

%Expression of analytical soln twilight zone method
syms t x y
% t = t;
% x = x;
% y = y;

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
  

% BC = fS.BC;

Nxg = fS.Nxg;
Nyg = fS.Nyg;
hx  = fS.hx;
hy  = fS.hy;
tic
Lp =  formLHSP(fS);
toc

fS.Lp = Lp;

% Lp2 = formLHSp1(fS);



if fS.directSolve == 0
    tic
    [Lpl,Lpu,Pp] = lu(Lp);
    fS.Lpl = Lpl;
    fS.Lpu = Lpu;
    fS.Pp  = Pp;
    toc
end

if  tExplicit == 0
    
    [Lu, Lv]   = formLHS(fS);
    fS.Lu = Lu;    
    fS.Lv = Lv;

    if fS.directSolve == 0
        
        [Lul,Luu,Up] = lu(Lu);
        [Lvl,Lvu,Vp] = lu(Lv);
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
[~,~,~,rhsuP1,rhsvP1] = VelocitySolverNew(tP1,tP1,count,fS);

fS.rhsuP1=rhsuP1;
fS.rhsvP1=rhsvP1;

if tOrder==4
    
    fS.UN = fS.UP2;
    fS.VN = fS.VP2;
    fS.PN = fS.PP2;
    fS.TemN = fS.TemP2;
    [~,~,~,rhsuP2,rhsvP2] = VelocitySolverNew(tP2,tP2,count,fS);
     
    fS.rhsuP2=rhsuP2;
    fS.rhsvP2=rhsvP2;

    fS.UN = fS.UP3;
    fS.VN = fS.VP3;
    fS.PN = fS.PP3;
    fS.TemN = fS.TemP3;
    [~,~,~,rhsuP3,rhsvP3] = VelocitySolverNew(tP3,tP3,count,fS);

    % check if the get rhs part is broken comparing to the previous version
    % okay 03/20/17
    %     [~,~,~,rhsuP3C,rhsvP3C] = VelocitySolver(tP3,tP3,count,fS);
    %     max(max(abs(rhsuP3C-rhsuP3)))
    %     max(max(abs(rhsvP3C-rhsvP3)))
    
    fS.rhsuP3=rhsuP3;
    fS.rhsvP3=rhsvP3;
    
end
%--------------------------------------------------------------------------
% Start time marching
%--------------------------------------------------------------------------
tic
%t1 = tC-dt;
% dtn = dt; 
% dte = dt;
n = 1;

% if fS.makeMovie == 1
%     surf(UN);
%     axis tight manual
%     set(gca,'nextplot','replacechildren'); 
% end

t2 = tC; 
t  = [textr,tP3,tP2,tP1,tC];

flagUpdate = 0;

while (t2 - tend) < eps-dt/4
    n = n+1;

    count = 1;
    t1 = t(end);
    t2 = t1 + dt;%the time of stage2 in each time step
    fS.t2= t2;

    t = [t t2];
    
    fS.tSpan = t(end-4:end);
    
    
    %% if dt is updated from the previous time step,
    %  recalculate coeffs and matrix..
    if flagUpdate == 1
        [fS.cA, fS.cBI, fS.cBE] = getBDFCoeff(fS);
        
        if  tExplicit == 0
            
            [Lu, Lv]   = formLHS(fS);
            fS.Lu = Lu;
            fS.Lv = Lv;
            
            if fS.directSolve == 0
                
                [Lul,Luu,Up] = lu(Lu);
                [Lvl,Lvu,Vp] = lu(Lv);
                fS.Lul = Lul;
                fS.Luu = Luu;
                fS.Up  = Up;
                
                fS.Lvl = Lvl;
                fS.Lvu = Lvu;
                fS.Vp  = Vp;
                
            end
            
        end
        if nUpdate < 4
            nUpdate = nUpdate+1;
        elseif nUpdate == 4
            flagUpdate = 0; % reset flagUpdate
        end
    end
     
    %% Predictor
    fS.UN = fS.UC;
    fS.VN = fS.VC;
    fS.PN = fS.PC;
    
    [count,U2,V2,rhsuC,rhsvC,maxgrad] = VelocitySolverNew(t1,t2,count,fS);
    
  Ut = fS.u(x,y,t2);
    Vt = fS.v(x,y,t2);
    Pt = fS.p(x,y,t2);
    Tt = fS.tem(x,y,t2);

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
    
    P2 = PressureSolverNew(t2,count,fS);
%     Pt = fS.p(x,y,t2);
%     max(max(abs(P2(3:end-2,3:end-2)-Pt(3:end-2,3:end-2))...
%          -abs(mean(mean(P2(3:end-2,3:end-2)-Pt(3:end-2,3:end-2))))))
%     
    fS.PN   = P2;
    
    for k = 1:correct
        
        [count,U2,V2] = VelocitySolverNew(t2,t2,count,fS);
%           U2 = fS.u(x,y,t2);
%     V2 = fS.v(x,y,t2);
    
        Tem2 = fS.tem(x,y,t2);
        
        fS.UN   = U2;
        fS.VN   = V2;
        fS.TemN = Tem2;
        
        P2 = PressureSolverNew(t2,count,fS);
        %P2 = fS.p(x,y,t2);
        
        fS.PN   = P2;
    end
    %

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
            
            if fS.tOrder==2
                
                stabilityReal = 1.6;
                stabilityImag = 1.28;
                
            elseif fS.tOrder==4
                
                if fS.tMethod == 1
                    
                    stabilityReal = 1.5;
                    stabilityImag = 1.16;
                    
                elseif fS.tMethod == 2
                    
                    stabilityReal = 1.;
                    stabilityImag = 1.2;
                    
                end
                
                
            end
            
            ilambda = max(max(U2))*(1+2/3)/hx + max(max(V2))*(1+2/3)/hy;
            if max(abs(ad41),abs(ad42)) > 0
                rlambda = 4*max(mu,al)*(1+1/3)/hx^2 + 4*max(mu,al)*(1+1/3)/hy^2 + 2*(ad41 + ad42*maxgrad)*16;
            else
                rlambda = 4*max(mu,al)*(1+1/3)/hx^2 + 4*max(mu,al)*(1+1/3)/hy^2;
            end
            dt    = sqrt(CFL/((ilambda/stabilityImag)^2 + (rlambda/stabilityReal)^2));
            
            %dt = dt *.5;
            tntemp = ceil((tend - t2)/dt);
            dt    = (tend - t2)/tntemp;
            fS.dt = dt;
            flagUpdate = 1;
            nUpdate    = 0;
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
            elseif fS.twilightZone == 9
                G  = f(5);
            end
            
            fx(x,y,t) = G + 0*x + 0*y + 0*t;
            fy(x,y,t) = 0*x + 0*y + + 0*t;
            
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
                fS.u  = matlabFunction(real(u));
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
        
        
        fx = simplify(dudt + u*dudx +v*dudy + dpdx - mu*dudx2 - mu*dudy2);
        fy = simplify(dvdt + u*dvdx +v*dvdy + dpdy - mu*dvdx2 - mu*dvdy2 - beta*g*(tem-tref));
        %pause
        dfxdx = diff(fx,x);
        dfydy = diff(fy,y);
        
        fxE = simplify(dudt + u*dudx +v*dudy + dpdx);
        fyE = simplify(dvdt + u*dvdx +v*dvdy + dpdy - beta*g*(tem-tref));
        
        fxI = simplify(- mu*dudx2 - mu*dudy2);
        fyI = simplify(- mu*dvdx2 - mu*dvdy2);
            
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
        fS.dpdy2  = matlabFunction(dpdy2);
        fS.dpdx2  = matlabFunction(dpdx2);
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

    function errorcalculation
        fprintf('at time %8.5f \n',t2)
        fprintf('l2 of U = %8.3e \n',err(1,n+1))
        fprintf('l2 of V = %8.3e \n',err(2,n+1))
        fprintf('l2 of T = %8.3e \n',err(3,n+1))
        fprintf('l2 of P = %8.3e \n',err(4,n+1))
        
    end

%if plotting == 1
    close all
    pl = figure('Visible','off');
    %pl = figure(1);
    %t = linspace(dt,tend,tn-1);
    label = 'Error vs t';
    
    t = t(4:end);
    
    plot(t,log10(err(1,:)),'LineWidth',2)
    hold on
    plot(t,log10(err(2,:)),'r--','LineWidth',2)
    %plot(t,log10(err(3,3:end)),'k','LineWidth',2)
    plot(t,log10(err(4,:)),'g','LineWidth',2)
    plot(t,log10(abs(diverr(:))),'m','LineWidth',2)

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
        view([90 90]);
        colorbar
        title('error V','FontSize',16);
        
        subplot(2,2,2);
        surf(x,y,U2-Ut)
        shading interp
        view([90 90]);
        colorbar
        title('error U','FontSize',16);
        %         subplot(2,2,3);
        %         surf(x,y,Tem2-Tt)
        %         title('error T','FontSize',16);
        subplot(2,2,3);
        
        surf(x,y,abs(P2-Pt) - abs(mean(mean(P2-Pt))))
        shading interp
        view([90 90]);
        colorbar
        title('error P','FontSize',16);
    end

pl=figure('visible','off');
debugPlot
print(pl, '-dpng',[fS.path 'Error' '.png']); 



if fS.makeMovie == 1
close(movieU);
close(movieV);
close(movieP);
end

%no temperature
for i = 1:2
    errf(i) = err(i,end);
end
errf(3)=err(4,end);

end



