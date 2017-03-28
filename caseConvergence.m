function caseConvergence(varargin)

% clearvars -except err rate
close all
nameappenSetup = '4th';
%8378 for test

cg3 = 0;

fS     =ins;
fS.tEnd=2.1;
fS.tOrder=4;
fS.tMethod=2;
fS.tExplicit=0;
fS.CFL=0;
fS.numberOfCorrector=1;
fS.numberOfAdapt=5000;
fS.extOrder=6;
fS.uvSwitch=0;
fS.CFLfix=.9;
fS.directSolve = 1;

fS.plotting = 0;
fS.makeMovie =0;

fS.BC = [1 1;1 1];
fS.twilightZone = 4; % 1 for trigSpolyT, 2 for poly2
if fS.twilightZone >= 7;
    fS.tw   = 0; % 1 for trigSpolyT, 2 for poly2
else
    fS.tw=1;
end
fS.cons = 1;         % 1 adv             2 divR         3 skew
fS.alpha = 0.01;
fS.WENO = 1;
fS.uw = 0;

gridx = [76 151];
gridy = [76 151];

gridx = [51 101];
gridy = [51 101];

%with cfl =.9, tw == 4, tEnd=1.1, gridx=gridy = 401; This takes roughly 9 hours with direct solve

fS.domain = [-1 1; -1 1];

if fS.twilightZone == 6
    fS.domain = [0 2*pi; 0 2*pi];
end

if fS.BC==2
    fS.al                 =  0.1;                                           % thermal diffusitity, alpha
    fS.K                  =  0;                                           % thermal conductivity K
    fS.mu                 =  0.1;                                           % Kinematic viscosity nu
    fS.tref               =  2;                                           % Treference T_{\infty};
    fS.beta               =  1;                                           % thermal expansion
    fS.g                  =  0;                                           % gravity g
    %par = [0,0,0,2,1,0];
else
    fS.al                 =  0.01;                                         % thermal diffusitity, alpha
    fS.K                  =  0;                                           % thermal conductivity K
    fS.mu                 =  0.01;                                         % Kinematic viscosity nu
    fS.tref               =  2;                                           % Treference T_{\infty};
    fS.beta               =  1;                                           % thermal expansion
    fS.g                  =  0;                                           % gravity g
    %par = [0.1,0,0.1,2,1,0];
end

    
%% get inline input argument

for i=1:nargin
    line = varargin{i};
    
    
    if(strncmp(line,'-Nx=',4))
        gridx = sscanf(varargin{i},'-Nx=%i');
    end

    if(strncmp(line,'-Ny=',4))
        gridy = sscanf(varargin{i},'-Ny=%i');
    end
    
    if(strncmp(line,'-tend=',6))
        fS.tEnd = sscanf(varargin{i},'-tend=%e');
    end
    if(strncmp(line,'-tOrder=',8))
        fS.tOrder = sscanf(varargin{i},'-tOrder=%i');
    end
    if(strncmp(line,'-fS.cons=',6))
        fS.cons = sscanf(varargin{i},'-fS.cons=%i');
    end
    if(strncmp(line,'-BC=',4))
      tempBC = sscanf(varargin{i},'-BC=%i');

      fS.BC  = ones(2,2)*tempBC;
    end
    if(strncmp(line,'-CFL=',5))
        fS.CFL = sscanf(varargin{i},'-CFL=%e');
    end
    if(strncmp(line,'-test=',6))
        fS.twilightZone = sscanf(varargin{i},'-test=%i');
    end
    if(strncmp(line,'-alpha=',7))
        fS.alpha = sscanf(varargin{i},'-alpha=%e');
    end
    if(strncmp(line,'-name',5))
        ls /Users/fanlongmeng/Documents/Research/INS4/tables
	nameappenSetup = input('appendix of the new file \n');
    end

end


if gridx(end) > 200 && gridy(end) > 200
    fS.directSolve = 1;    
end


if fS.cons == 1
    for i=1:nargin
        line = varargin{i};
        if (strncmp(line,'-WENO=',6))
            fS.WENO = sscanf(varargin{i},'-WENO=%i');
        end
        if(strncmp(line,'-uw=',4))
            fS.uw = sscanf(varargin{i},'-uw=%i');
        end
    end
end

%testfunction
syms t x y

switch fS.twilightZone  % u::1 v::2 p::3 tem::4
    case 1
        f(1) = (sin(pi*x))^2*sin(2*pi*y)*(1+t);
        f(2) = -sin(2*pi*x)*(sin(pi*y))^2*(1+t);
        f(3) = (sin(pi*x)*sin(pi*y))*(1+t);
        f(4) = (sin(pi*y)*sin(pi*x))*(1+t);
    case 2
        if fS.tOrder==2
            f(1) = (x^2+2*x*y+y^2)*(1+2*t+2*t^2);
            f(2) = (-x^2-2*x*y-y^2)*(1+2*t+2*t^2);
            f(3) = (x^2+1/2*x*y+y^2-1)*(1+2*t+2*t^2);
            f(4) = (x^2+2*x*y+y^2)*(1+2*t+2*t^2);
        elseif fS.tOrder==4
            if fS.tMethod == 1
                
                            f(1) = (x^2+2*x*y+y^2)*(1+2*t+2*t^2+1/3*t^3);
            f(2) = (-x^2-2*x*y-y^2)*(1+2*t+2*t^2+1/3*t^3);
            f(3) = (x^2+1/2*x*y+y^2-1)*(1+2*t+2*t^2+1/3*t^3);
            f(4) = (x^2+2*x*y+y^2)*(1+2*t+2*t^2+1/3*t^3);
            
            elseif fS.tMethod==2
                
                %             f(1) = (x^2+2*x*y+y^2)*(1+2*t+2*t^2+1/3*t^3+1/4*t^4);
                %             f(2) = (-x^2-2*x*y-y^2)*(1+2*t+2*t^2+1/3*t^3+1/4*t^4);
                %             f(3) = (x^2+1/2*x*y+y^2-1)*(1+2*t+2*t^2+1/3*t^3+1/4*t^4);
                %             f(4) = (x^2+2*x*y+y^2)*(1+2*t+2*t^2+1/3*t^3+1/4*t^4);
                
                f(1) = (x^2+2*x*y+y^2)*(1+2*t+2*t^2);
                f(2) = (-x^2-2*x*y-y^2)*(1+2*t+2*t^2);
                f(3) = (x^2+1/2*x*y+y^2-1)*(1+2*t+2*t^2);
                f(4) = (x^2+2*x*y+y^2)*(1+2*t+2*t^2);
                
%                 f(1) = ( x + y)*(1) ;
%                 f(2) = (-x - y)*(1) ;
%                 f(3) = (x+y)*(1)*eps + 1;
%                 f(4) = (x+y)*(1);
                
            end
        end
    case 3
        f(1) = (x^3+y^3-3*x*y^2)*(1+t);
        f(2) = (x^3+y^3-3*x^2*y)*(1+t);
        f(3) = (x^3+1/2*x*y+y^3-1)*(1+t);
        f(4) = (2+x+y/2+x^2/2+y^2/4)*(1+t);
    case 4 %% trig
        f(1) = (sin(pi*x))^2*sin(2*pi*y)*(cos(pi*t*0.2));
        f(2) = -sin(2*pi*x)*(sin(pi*y))^2*(cos(pi*t*0.2));
        f(3) = (sin(pi*x)*sin(pi*y))*(cos(pi*t*0.2));
        f(4) = eps.*(sin(pi*y)*sin(pi*x))*(cos(pi*t*0.2)+1);
    case 5  %% trig u,v always larger to 0
        f(1) = (sin(pi*x))^2*sin(2*pi*y)*(cos(pi*t)+2)+4;
        f(2) = -sin(2*pi*x)*(sin(pi*y))^2*(cos(pi*t)+2)+4;
        f(3) = (sin(pi*x)*sin(pi*y)+2)*(cos(pi*t)+2)+4;
        f(4) = eps.*(sin(pi*y)*sin(pi*x))*(cos(pi*t)+2)+4;
    case 6  %% taylor green in domain [0,2*pi] with periodicity 
         fS.beta               =  0;
        f(1) =  sin(x)*cos(y)*(exp(-2*fS.mu*t));
        f(2) =  sin(y)*cos(x)*(-exp(-2*fS.mu*t));
        f(3) = 1/4*(cos(2*x) + cos(2*y))*(exp(-4*fS.mu*t));
        f(4) = eps.*(sin(pi*y)*sin(pi*x))*(cos(pi*t)+2)+2;
    case 7
        G = 1;
        H = fS.domain(2,2);
        mu = fS.mu;
        f(1) = eps*x + eps*y + eps*t+ G*H^2/(2*mu)*(1 - y.^2/H^2);
        f(2) = eps*x + eps*y + eps*t;
        f(3) = eps*x + eps*y + eps*t;
        f(4) = eps*x + eps*y + eps*t;
    case 8 
        G = 1;
        H = fS.domain(2,2);
        delta = 0.25;
        
        f(1) = eps*x + eps*y + eps*t+ ((y>0).*G*H^2/(2)*(1 - y.^2/H^2)) ...
            + ((y<=0).*(-G)*H^2/(2)*(1 - y.^2/H^2));
        f(2) = delta*(sin(-pi*x)).*exp(-10*(y.^2/H^2))  + eps*t;
        f(3) = eps*x + eps*y + eps*t;
        f(4) = eps*x + eps*y + eps*t;    
end

for i = 1:length(gridx)
    
    fS.Nx = gridx(i);
    fS.Ny = gridy(i);
    h(i) = (fS.domain(1,2)-fS.domain(1,1))/(fS.Nx-1);

    
    %tvar = [CFL,tend,correct,nadapt,tOrder,fourthOrderT];
    
    fS.ad41 = 0;
    fS.ad42 = 0;
    fS.ad21 = 0;
    fS.ad22 = 0;
    
    %stabilization = [ad41,ad42,ad21,ad22,WENO,uw];
    
    nameappen=['T',num2str(fS.tOrder)];
    
    if fS.tOrder==4
        nameappen=[nameappen,'FT',num2str(fS.tMethod)]; %FourthTime
    end
    
    if (fS.mu == 0 && fS.al == 0)
        nameappen=[nameappen,'Nu',num2str(0),'S',nameappenSetup];
    end
    
    nameappen=[nameappen,'S',nameappenSetup];

            
    fS.path = ['/Users/fanlongmeng/Documents/Research/INS4/figures/figure' nameappen];
    if cg3==1
           fS.path = ['/home/mengf5/fins/figures/figure' nameappen]; 
    end
        
    
    err(i,:) = main(f,fS);
    
end

clc

for i = 1:length(gridx)
    fprintf('mesh size = %8.5f err = %e %e %e \n',h(i),err(i,:))
end

for j = 1:length(gridx)-1
    rate(j,1:3) = log(err(j,:)./err(j+1,:))./log(h(j)/h(j+1));
end

%end
time = char(datestr(now,'mmmm dd, yyyy HH:MM:SS AM'));
%save tableContent.mat err rate

pathOutput = ['/Users/fanlongmeng/Documents/Research/INS4/setups/setup' nameappen '.tex'];
if cg3==1

pathOutput = ['/home/mengf5/fins/setups/setup' nameappen '.tex'];
end

FID = fopen(pathOutput, 'w');
fprintf(FID, 'date saved at %s \\\\ \n',time);
fprintf(FID, 'x-domain is $\\in [%3.2f, %3.2f]$ \\\\ \n',fS.domain(1,:));
fprintf(FID, 'y-domain is $\\in [%3.2f, %3.2f]$ \\\\ \n',fS.domain(2,:));
fprintf(FID, 'gridx = ');
fprintf(FID,[repmat(' %d ', 1, length(gridx)) '\\\\ \n'], gridx');
fprintf(FID, 'gridy = ');
fprintf(FID,[repmat(' %d ', 1, length(gridx)) '\\\\ \n'], gridx');

if fS.CFL == 0
    fprintf(FID, 'The simulation was calculated witn fixed dt based on the initial condition and CFL = %3.2f \\\\ \n',fS.CFLfix);
else
    fprintf(FID, 'CFL = %i \\\\ \n',fS.CFL);
    fprintf(FID, 'dt is updated at every %i th step\\\\ \n',fS.numberOfAdapt);
end
fprintf(FID, 'final time = %3.2f \\\\ \n',fS.tEnd);
fprintf(FID, 'there is %i correcter step used in time integration \\\\ \n',fS.numberOfCorrector);

fprintf(FID, 'thermal diffusitity $\\alpha = %3.2f$ \\\\ \n',fS.al);
fprintf(FID, 'kinematic viscosity $\\nu= %3.2f$ \\\\ \n',fS.mu);
fprintf(FID, 'reference temperature $T_{\\infty} = %3.2f $ \\\\ \n',fS.tref);
fprintf(FID, 'thermal expansion rate $\\beta = %3.2f $ \\\\ \n',fS.beta);
fprintf(FID, 'gravity $g = %3.2f $ \\\\ \n',fS.g);

if fS.BC == 1
    fprintf(FID, 'Dirichlet Boundary condition is imposed at the boundary \\\\ \n');
elseif fS.BC == 2
    fprintf(FID, 'Periodic Boundary condition is imposed at the boundary \\\\ \n');
elseif fS.BC == 3
    fprintf(FID, 'Exact Boundary condition is imposed at the boundary and ghost points \\\\ \n');
elseif fS.BC == 6
    fprintf(FID, 'periodic boundary condition on x-direction, Dirichlet on y-direction \\\\ \n');

end

fprintf(FID, 'Test function is: \\\\ \n');
fprintf(FID, '$u = %s$ \\\\ \n',char(f(1)));
fprintf(FID, '$v = %s$ \\\\ \n',char(f(2)));
fprintf(FID, '$p = %s$ \\\\ \n',char(f(3)));
fprintf(FID, '$tem = %s$ \\\\ \n',char(f(4)));

if abs(fS.ad41) > 0 || abs(fS.ad42) > 0 || abs(fS.ad21) > 0 || abs(fS.ad22) > 0
    fprintf(FID, 'artificial dissipation method is used with: \\\\ \n');
    fprintf(FID, '$ad41 == %3.2f$, $ad42 == %3.2f$, $ad21 == %3.2f$, $ad22 == %3.2f$ \\\\ \n', fS.ad41,fS.ad42,fS.ad21,fS.ad22);
end

if fS.WENO > 0
    if fS.uw > 0
        fprintf(FID, 'Upwind biased scheme is used \\\\ \n');
    else
        fprintf(FID, 'WENO scheme is used \\\\ \n');
    end
end

fprintf(FID, 'The constant for divergence damping is $\\alpha_d == %3.2f$ \\\\ \n',fS.alpha);

switch fS.cons
    case 1
        fprintf(FID, 'the convective term is discretized with adv \\\\ \n');
    case 2
        fprintf(FID, 'the convective term is discretized with divRs \\\\ \n');
    case 3
        fprintf(FID, 'the convective term is discretized with skew \\\\ \n');
end

fclose(FID);


pathTable = ['/Users/fanlongmeng/Documents/Research/INS4/tables/table' nameappen '.tex'];

if cg3==1

pathTable = ['/home/mengf5/fins/tables/table' nameappen '.tex'];

end

FID = fopen(pathTable, 'w');
fprintf(FID, '\\begin{centering}  \n');
%fprintf(FID, '\\begin{tabular}{c|l|l|l|c} \n');
fprintf(FID, '\\begin{tabular}{c|l|l|c} \n');
fprintf(FID, '\\hline \n');
fprintf(FID, '\\textrm{} &  $|U-U_{exact}|$ & $|V-V_{exact}|$ &  $|P-P_{exact}|$  \\\\ \n');
%fprintf(FID, '\\textrm{} &  $|U-U_{exact}|$ & $|V-V_{exact}|$ & $|T-T_{exact}|$ & $|P-P_{exact}|$  \\\\ \n');
fprintf(FID, '\\hline \n');

for k=1:length(err(:,1))
    fprintf(FID, '$h = %5.4f$ & %8.2e & %8.2e &  %8.2e  \\\\ ', h(k),err(k,:));
    if k==length(err(:,1))
        fprintf(FID, '\\hline ');
    end
    fprintf(FID, '\n');
end
if length(gridx) ~= 1
   for k = 1:length(gridx)-1
       fprintf(FID, '$rate$ & %3.2f & %3.2f &  %3.2f  \\\\ \n ',rate(k,:));
   end
       fprintf(FID, '\\hline \n');
end

fprintf(FID, '\\end{tabular}\n');
fprintf(FID, '\\end{centering} \n');
fclose(FID);
