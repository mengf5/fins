classdef ins < handle
    %% the ins class include properties w.r.t
    % time integration, spatial, scheme coefficients,
    % and physics coefficients.
    % the functions def here are:
    %     setup(fS)
    %     getDt(fS,u,v)
   properties     
      %% time integration related variables 
      tExplicit          =  1;
      tEnd               =  1;
      tOrder             =  2;
      tMethod            =  1;
      fLCBDF             =  0;
      cA                 = []; % coefficients of BDF of previous u(LHS)
      cBI                = []; % coefficients of BDF of fI(RHS)
      cBE                = []; % coefficients of BDF of fE(RHS)      
      CFL                =  0;
      CFLfix             =  .9;
      numberOfCorrector  =  0;
      numberOfAdapt      =  0;
%      preTime            = [];
%      corrTime           = [];
%      imTime             = [];
      dt; dtn; dte;
      tSpan              = [];
      t2; % this is the current time 
      
      %% LHS matrix for U,V,P
      Lul; Luu; Up;
      Lvl; Lvu; Vp;
      Lpl; Lpu; Pp;
      Lp; Lu; Lv;
      directSolve  = 0;
      
      %% spatial related variables
      Nx                 = 10;
      Ny                 = 10;
      BC                 = [2 2;2 2]; % BC(axis,side)  
                                      % 1. noslip ;    2. periodic ;
                                      % 3. dirichlet ; 4. inlet; 
                                      % 5. outlet;
                                      % 6. Dirichlet
      domain             = [-1 1; -1 1];
      Nxg;      
      Nyg;
      hx;
      hy;
      Kbx;Kby;albx;alby;mubx;muby;
      bcp=[0 0;0 0];
      x;y;   
      gL         = 2;
      ia;ib;ja;jb;
      
      %% scheme related variables ( see cginRef for more info.)
      twilightZone       =  1;                                              % choose which twilightZone forcing function
      tw                 =  1;
      alpha              =  1;                                              % constant coeffcient of divergence damping
      cons               =  1;                                              % if using conservative variables for advection terms, 1 is the only working scheme!
      ad41=0; ad42=0; ad21=0; ad22=0;                                       % coefficients for strong artifitial dissipation
      WENO               =  1;                                              % if using BWENO scheme 
      uw                 =  0;                                              % if using thrid order upwinding scheme (as a milestone scheme for stability testing)
      extOrder           =  6;                                              % order of extrapolation on the ghost points
      uvSwitch           =  1;                                              % switch the order of pre-corr for solving u and v
      %% physics coefficients
      al                 =  0;                                              % thermal diffusitity, alpha
      K                  =  0;                                              % thermal conductivity K
      mu                 =  0;                                              % Kinematic viscosity nu
      tref               =  0;                                              % Treference T_{\infty}; 
      beta               =  0;                                              % thermal expansion
      g                  =  9.81;                                           % gravity g
      
      %% supporting I/O 
      path;
      plotting           = 0;
      makeMovie          = 0;
      
      %% solutions 
      UN;   UC;   UP1;   UP2;   UP3;
      VN;   VC;   VP1;   VP2;   VP3;
      PN;   PC;   PP1;   PP2;   PP3;    Pextr;
      TemN; TemC; TemP1; TemP2; TemP3;
      
      rhsuP1; rhsuP2; rhsuP3;
      rhsvP1; rhsvP2; rhsvP3;
      rhstP1; rhstP2; rhstP3;
      
      %% function handles for forcing and bcs
      u;v;p;tem;
      fx;fy;ftem; dfxdx;dfydy;
      dudt;dvdt; dubcdt;dvbcdt;
      dudx; dvdx; dudy; dvdy;
      dvdy2; dvdy3; dvdx2; dudx2; dudy2;
      fxE;fyE;fxI;fyI;
      dpdx;dpdy;dudxy;dvdxy;dvdty;dvdx2y;
      
      
   end
   methods
       function setup(fS)
           
           
           fS.Nxg = fS.Nx + 4; % # of x-grid + ghost points
           fS.Nyg = fS.Ny + 4; % # of y-grid + ghost points
           fS.ia  = fS.gL+1;
           fS.ib  = fS.Nxg - fS.gL;
           fS.ja  = fS.gL+1;
           fS.jb  = fS.Nyg - fS.gL;
      
           %geometry
           xs = fS.domain(1,1); % starting point in x
           xe = fS.domain(1,2); % end point in x
           
           ys = fS.domain(2,1); % starting point in y
           ye = fS.domain(2,2); % end point in y
           
           fS.hx  = (xe-xs)/(fS.Nx-1);
           fS.hy  = (ye-ys)/(fS.Ny-1);
           
           fS.Kbx  = fS.K/(12*fS.hx);
           fS.Kby  = fS.K/(12*fS.hy);
           
           fS.albx = fS.al/(12*fS.hx^2);
           fS.alby = fS.al/(12*fS.hy^2);
           
           fS.mubx = fS.mu/(12*fS.hx^2);
           fS.muby = fS.mu/(12*fS.hy^2);
           
           fS.bcp(1,1)= 3;
           fS.bcp(1,2)= fS.Nxg-2;
           fS.bcp(2,1)= 3;
           fS.bcp(2,2)= fS.Nyg-2;
           
           [x,y] = meshgrid(linspace(xs-2*fS.hx,xe+2*fS.hx,fS.Nxg),linspace(ys-2*fS.hy,ye+2*fS.hy,fS.Nyg));
           fS.x = x';
           fS.y = y';
           
       end
       
       function getDt(fS)
           
           %r = [insVar.CFL]+[insVar.Nx];
           
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
           
           % get the initial U0 and V0
           if fS.twilightZone == 7
               G = 1;
               H = fS.domain(2,2);
               U0 = G*H^2/(2*fS.mu)*(1 - fS.y.^2/H^2);
           else
               U0   = fS.u(fS.x,fS.y,0);
           end
           V0   = fS.v(fS.x,fS.y,0);
           
           ilambda = max(max(U0))*(1+2/3)/fS.hx + max(max(V0))*(1+2/3)/fS.hy;
           
           if fS.tExplicit==0
               rlambda = 0;
               
           elseif fS.tExplicit==1
               
               if max(abs(fS.ad41),abs(fS.ad42)) > 0
                   gradu = gradient(U0);
                   gradv = gradient(V0);
                   maxgrad = max(max(max(gradu,gradv)));
                   rlambda = 4*max(fS.mu,fS.al)*(1+1/3)/fS.hx^2 + 4*max(fS.mu,fS.al)*(1+1/3)/fS.hy^2 + 2*(fS.ad41 + fS.ad42*maxgrad)*16;
               else
                   rlambda = 4*max(fS.mu,fS.al)*(1+1/3)/fS.hx^2 + 4*max(fS.mu,fS.al)*(1+1/3)/fS.hy^2;
               end
           end
           
           %nadapt = 10;
           if fS.CFL ~=0
               tempCFL = fS.CFL;
           else
               tempCFL = fS.CFLfix;
           end
           
           fS.dt = sqrt(tempCFL/((ilambda/stabilityImag)^2 + (rlambda/stabilityReal)^2));
           
           if fS.tMethod == 1
               stabilityRealViscous = 3;
               rlambda = 4*max(fS.mu,fS.al)*(1+1/3)/fS.hx^2 + 4*max(fS.mu,fS.al)*(1+1/3)/fS.hy^2;
               dtViscous = sqrt(fS.CFLfix/( (rlambda/stabilityRealViscous)^2));
               if dtViscous < fS.dt
                   fS.dt = dtViscous;
               end
           end
           %t = dt;
           tn = ceil(fS.tEnd/fS.dt);       % number of time-step
           fS.dt = fS.tEnd/tn;
           
           
           if fS.tOrder==2
               
%               fS.preTime  = [3/2*fS.dt,-1/2*fS.dt];
%               fS.corrTime = [1/2*fS.dt, 1/2*fS.dt];
           
           elseif fS.tOrder==4
               
               if fS.tMethod==1
                   %fS.preTime =
                   %[55/24*fS.dt,-59/24*fS.dt,37/24*fS.dt,-3/8*fS.dt]; this
                   % is ab 4
                   %fS.preTime = [23/12*fS.dt,-4/3*fS.dt,5/12*fS.dt,0]; 
                   % this is ab3
                   
                   %fS.corrTime =[3/8*fS.dt,19/24*fS.dt,-5/24*fS.dt,1/24*fS.dt];
                   
               elseif fS.tMethod==2
                   
                   fS.cBE = 12/25*[4,-6,4,-1]*fS.dt;
                   %fS.preTime = fS.cBE;
                   
                   fS.cBI = 12/25*fS.dt;
                   %fS.corrTime = [fS.cBI,0,0,0];
                   
                   fS.cA  = 12/25*[-4,3,-4/3,1/4];
                   
               end
           end
           
           if fS.tExplicit==0
               if fS.tOrder==2
                   %fS.imTime = [1/2*fS.dt, 1/2*fS.dt];
               elseif fS.tOrder==4
                   if fS.tMethod==1                       
                       %fS.imTime = [3/8*fS.dt,19/24*fS.dt,-5/24*fS.dt,1/24*fS.dt];
                   elseif fS.tMethod==2
                       fS.cBI = 12/25*fS.dt;
                       %fS.imTime = fS.cBI;
                   end
               end
           end
           
           
           fS.dtn = fS.dt; 
           fS.dte = fS.dt;
     
       end
       %        function r = multiplyBy(obj,n)
       %            r = [obj.Ny] * n;
       %        end
   

  end 
end