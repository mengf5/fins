function imexTimeStepping
close all
xv = linspace(-.1,.1,201);
yv = linspace(-1,1,201);

[X,Y] = meshgrid(xv,yv); %generate grids for contour
Z = X+1i*Y; % complex number

order = 4;
fixStep=0;
localOrder = 5;
addTerm = 'R';
newTimeLevel = 2;
fLCBDF =0 ;

%% IMEXBDF

for k = 1:5
    
    r = 1 + (k-3)*0.2;
    %r = 1;
    
    lC = ['c','b','r','m','k'];
    
    if order == 2
        
        tnp1 = 0.1;
        tn   = 0.09;
        tnm1 = 0.075;
        
        r1  = tn/tnp1;
        r2  = tnm1/tnp1;
        fLC = 2/3;
        a1  =  -(2*r1 - 3*r2 + 1)/(3*(r1 - r2));
        a2  = -(r1 - 1)/(3*(r1 - r2));
        
        b0  = fLC*(1-r1);
        b1  = (2*r1 - 5*r2 + 3)/(6*(r1 - r2));
        b2  = (2*r1 + r2 - 3)/(6*(r1 - r2));
        
    elseif order==4
        
        dt  = 0.01;
        dtN = dt*r;
        %r  = 1;
        
        tnm4 = .1;
        tnm3 = tnm4 + dt;
        
        if newTimeLevel == 4
            dt = dtN;
        end
        tnm2 = tnm3 + dt;
        
        if newTimeLevel == 3
            dt = dtN;
        end
        tnm1 = tnm2 + dt;
        
        if newTimeLevel == 2
            dt = dtN;
        end
        tn   = tnm1 + dt;
        
        if newTimeLevel == 1
            dt = dtN;
        end
        tnp1 = tn   + dt;
        
        r1 = tn/tnp1;
        r2 = tnm1/tnp1;
        r3 = tnm2/tnp1;
        r4 = tnm3/tnp1;
        r5 = tnm4/tnp1;
        
        fLC = 12/25;
        
        if localOrder == 4
            
            a1=    (r2 - 36*r1 + r3 + r4 + 24*r1*r2 + 24*r1*r3 + 24*r1*r4 - 13*r2*r3 - 13*r2*r4 - 13*r3*r4 - 12*r1*r2*r3 - 12*r1*r2*r4 - 12*r1*r3*r4 + 25*r2*r3*r4 + 11)/(25*(r1 - r2)*(r1 - r3)*(r1 - r4));
            a2=                                                -(r3 - 35*r1 + r4 + 11*r1*r3 + 11*r1*r4 - 13*r3*r4 - 12*r1^2*r3 - 12*r1^2*r4 + 24*r1^2 + 13*r1*r3*r4 + 11)/(25*(r1 - r2)*(r2 - r3)*(r2 - r4));
            a3=                                      (r2 - 35*r1 + r4 + 11*r1*r2 + 11*r1*r4 - 13*r2*r4 - 12*r1^2*r2 - 12*r1^2*r4 + 24*r1^2 + 13*r1*r2*r4 + 11)/(25*(r3 - r4)*(r1*r2 - r1*r3 - r2*r3 + r3^2));
            a4 =-(r2 - 35*r1 + r3 + 11*r1*r2 + 11*r1*r3 - 13*r2*r3 - 12*r1^2*r2 - 12*r1^2*r3 + 24*r1^2 + 13*r1*r2*r3 + 11)/(25*(r1*r4^2 + r2*r4^2 + r3*r4^2 - r4^3 + r1*r2*r3 - r1*r2*r4 - r1*r3*r4 - r2*r3*r4));
            a5 = 0;
            
            
            b0 = fLC*(1-r1);
            b1=                                           (36*r1 - 37*r2 - 37*r3 - 37*r4 - 24*r1*r2 - 24*r1*r3 - 24*r1*r4 + 49*r2*r3 + 49*r2*r4 + 49*r3*r4 + 12*r1*r2*r3 + 12*r1*r2*r4 + 12*r1*r3*r4 - 61*r2*r3*r4 + 25)/(100*(r1 - r2)*(r1 - r3)*(r1 - r4));
            b2=                                                -(11*r2 - 12*r1 - 37*r3 - 37*r4 - 24*r1*r2 + 24*r1*r3 + 24*r1*r4 + r2*r3 + r2*r4 + 49*r3*r4 + 12*r1*r2*r3 + 12*r1*r2*r4 - 36*r1*r3*r4 - 13*r2*r3*r4 + 25)/(100*(r1 - r2)*(r2 - r3)*(r2 - r4));
            b3=                                      (11*r3 - 37*r2 - 12*r1 - 37*r4 + 24*r1*r2 - 24*r1*r3 + 24*r1*r4 + r2*r3 + 49*r2*r4 + r3*r4 + 12*r1*r2*r3 - 36*r1*r2*r4 + 12*r1*r3*r4 - 13*r2*r3*r4 + 25)/(100*(r3 - r4)*(r1*r2 - r1*r3 - r2*r3 + r3^2));
            b4 = -(11*r4 - 37*r2 - 37*r3 - 12*r1 + 24*r1*r2 + 24*r1*r3 - 24*r1*r4 + 49*r2*r3 + r2*r4 + r3*r4 - 36*r1*r2*r3 + 12*r1*r2*r4 + 12*r1*r3*r4 - 13*r2*r3*r4 + 25)/(100*(r1*r4^2 + r2*r4^2 + r3*r4^2 - r4^3 + r1*r2*r3 - r1*r2*r4 - r1*r3*r4 - r2*r3*r4));
            
        end
        
        if localOrder == 5
            
            if fLCBDF == 0
                
                
                a1 =                                                    (r2^2*r3^2*r4^2 - 2*r2^2*r3^2*r4 + r2^2*r3^2 - 2*r2^2*r3*r4^2 + 4*r2^2*r3*r4 - 2*r2^2*r3 + r2^2*r4^2 - 2*r2^2*r4 + r2^2 - 2*r2*r3^2*r4^2 + 4*r2*r3^2*r4 - 2*r2*r3^2 + 4*r2*r3*r4^2 - 8*r2*r3*r4 + 4*r2*r3 - 2*r2*r4^2 + 4*r2*r4 - 2*r2 + r3^2*r4^2 - 2*r3^2*r4 + r3^2 - 2*r3*r4^2 + 4*r3*r4 - 2*r3 + r4^2 - 2*r4 + 1)/((r1 - r2)*(r1 - r3)*(r1 - r4)*(3*r1 + 3*r2 + 3*r3 + 3*r4 - 2*r1*r2 - 2*r1*r3 - 2*r1*r4 - 2*r2*r3 - 2*r2*r4 - 2*r3*r4 + r1*r2*r3 + r1*r2*r4 + r1*r3*r4 + r2*r3*r4 - 4));
                a2 =                                                   -(r1^2*r3^2*r4^2 - 2*r1^2*r3^2*r4 + r1^2*r3^2 - 2*r1^2*r3*r4^2 + 4*r1^2*r3*r4 - 2*r1^2*r3 + r1^2*r4^2 - 2*r1^2*r4 + r1^2 - 2*r1*r3^2*r4^2 + 4*r1*r3^2*r4 - 2*r1*r3^2 + 4*r1*r3*r4^2 - 8*r1*r3*r4 + 4*r1*r3 - 2*r1*r4^2 + 4*r1*r4 - 2*r1 + r3^2*r4^2 - 2*r3^2*r4 + r3^2 - 2*r3*r4^2 + 4*r3*r4 - 2*r3 + r4^2 - 2*r4 + 1)/((r1 - r2)*(r2 - r3)*(r2 - r4)*(3*r1 + 3*r2 + 3*r3 + 3*r4 - 2*r1*r2 - 2*r1*r3 - 2*r1*r4 - 2*r2*r3 - 2*r2*r4 - 2*r3*r4 + r1*r2*r3 + r1*r2*r4 + r1*r3*r4 + r2*r3*r4 - 4));
                a3 =                                         (r1^2*r2^2*r4^2 - 2*r1^2*r2^2*r4 + r1^2*r2^2 - 2*r1^2*r2*r4^2 + 4*r1^2*r2*r4 - 2*r1^2*r2 + r1^2*r4^2 - 2*r1^2*r4 + r1^2 - 2*r1*r2^2*r4^2 + 4*r1*r2^2*r4 - 2*r1*r2^2 + 4*r1*r2*r4^2 - 8*r1*r2*r4 + 4*r1*r2 - 2*r1*r4^2 + 4*r1*r4 - 2*r1 + r2^2*r4^2 - 2*r2^2*r4 + r2^2 - 2*r2*r4^2 + 4*r2*r4 - 2*r2 + r4^2 - 2*r4 + 1)/((r3 - r4)*(r1*r2 - r1*r3 - r2*r3 + r3^2)*(3*r1 + 3*r2 + 3*r3 + 3*r4 - 2*r1*r2 - 2*r1*r3 - 2*r1*r4 - 2*r2*r3 - 2*r2*r4 - 2*r3*r4 + r1*r2*r3 + r1*r2*r4 + r1*r3*r4 + r2*r3*r4 - 4));
                a4 =-(r1^2*r2^2*r3^2 - 2*r1^2*r2^2*r3 + r1^2*r2^2 - 2*r1^2*r2*r3^2 + 4*r1^2*r2*r3 - 2*r1^2*r2 + r1^2*r3^2 - 2*r1^2*r3 + r1^2 - 2*r1*r2^2*r3^2 + 4*r1*r2^2*r3 - 2*r1*r2^2 + 4*r1*r2*r3^2 - 8*r1*r2*r3 + 4*r1*r2 - 2*r1*r3^2 + 4*r1*r3 - 2*r1 + r2^2*r3^2 - 2*r2^2*r3 + r2^2 - 2*r2*r3^2 + 4*r2*r3 - 2*r2 + r3^2 - 2*r3 + 1)/((r1*r4^2 + r2*r4^2 + r3*r4^2 - r4^3 + r1*r2*r3 - r1*r2*r4 - r1*r3*r4 - r2*r3*r4)*(3*r1 + 3*r2 + 3*r3 + 3*r4 - 2*r1*r2 - 2*r1*r3 - 2*r1*r4 - 2*r2*r3 - 2*r2*r4 - 2*r3*r4 + r1*r2*r3 + r1*r2*r4 + r1*r3*r4 + r2*r3*r4 - 4));
                
                b0 =                                                                                                                                                                                                                                                                                                                                                    (r2 + r3 + r4 - r2*r3 - r2*r4 - r3*r4 + r2*r3*r4 - 1)/(3*r1 + 3*r2 + 3*r3 + 3*r4 - 2*r1*r2 - 2*r1*r3 - 2*r1*r4 - 2*r2*r3 - 2*r2*r4 - 2*r3*r4 + r1*r2*r3 + r1*r2*r4 + r1*r3*r4 + r2*r3*r4 - 4);
                
                
                b1 =                                                                                                                                            -(r2^2*r3^2*r4^2 - 2*r2^2*r3^2*r4 + r2^2*r3^2 - 2*r2^2*r3*r4^2 + 4*r2^2*r3*r4 - 2*r2^2*r3 + r2^2*r4^2 - 2*r2^2*r4 + r2^2 - 2*r2*r3^2*r4^2 + 4*r2*r3^2*r4 - 2*r2*r3^2 + 4*r2*r3*r4^2 - 8*r2*r3*r4 + 4*r2*r3 - 2*r2*r4^2 + 4*r2*r4 - 2*r2 + r3^2*r4^2 - 2*r3^2*r4 + r3^2 - 2*r3*r4^2 + 4*r3*r4 - 2*r3 + r4^2 - 2*r4 + 1)/((r1 - r2)*(r1 - r3)*(r1 - r4)*(3*r1 + 3*r2 + 3*r3 + 3*r4 - 2*r1*r2 - 2*r1*r3 - 2*r1*r4 - 2*r2*r3 - 2*r2*r4 - 2*r3*r4 + r1*r2*r3 + r1*r2*r4 + r1*r3*r4 + r2*r3*r4 - 4));
                b2 =                                    -(r1 + r2 + 2*r3 + 2*r4 - r3^2*r4^2 - r1*r2 - 2*r1*r3 - 2*r1*r4 - 2*r2*r3 - 2*r2*r4 - 4*r3*r4 + r1*r3^2 + r1*r4^2 + r2*r3^2 + r2*r4^2 + 2*r3*r4^2 + 2*r3^2*r4 - r3^2 - r4^2 - r1*r2*r3^2 - r1*r2*r4^2 - 2*r1*r3*r4^2 - 2*r1*r3^2*r4 - 2*r2*r3*r4^2 - 2*r2*r3^2*r4 + r1*r3^2*r4^2 + r2*r3^2*r4^2 + 2*r1*r2*r3 + 2*r1*r2*r4 + 4*r1*r3*r4 + 4*r2*r3*r4 - 4*r1*r2*r3*r4 + 2*r1*r2*r3*r4^2 + 2*r1*r2*r3^2*r4 - r1*r2*r3^2*r4^2 - 1)/((r1 - r2)*(r2 - r3)*(r2 - r4)*(3*r1 + 3*r2 + 3*r3 + 3*r4 - 2*r1*r2 - 2*r1*r3 - 2*r1*r4 - 2*r2*r3 - 2*r2*r4 - 2*r3*r4 + r1*r2*r3 + r1*r2*r4 + r1*r3*r4 + r2*r3*r4 - 4));
                b3 =                          (r1 + 2*r2 + r3 + 2*r4 - r2^2*r4^2 - 2*r1*r2 - r1*r3 - 2*r1*r4 - 2*r2*r3 - 4*r2*r4 - 2*r3*r4 + r1*r2^2 + r1*r4^2 + r2^2*r3 + 2*r2*r4^2 + 2*r2^2*r4 + r3*r4^2 - r2^2 - r4^2 - r1*r2^2*r3 - 2*r1*r2*r4^2 - 2*r1*r2^2*r4 - r1*r3*r4^2 - 2*r2*r3*r4^2 - 2*r2^2*r3*r4 + r1*r2^2*r4^2 + r2^2*r3*r4^2 + 2*r1*r2*r3 + 4*r1*r2*r4 + 2*r1*r3*r4 + 4*r2*r3*r4 - 4*r1*r2*r3*r4 + 2*r1*r2*r3*r4^2 + 2*r1*r2^2*r3*r4 - r1*r2^2*r3*r4^2 - 1)/((r3 - r4)*(r1*r2 - r1*r3 - r2*r3 + r3^2)*(3*r1 + 3*r2 + 3*r3 + 3*r4 - 2*r1*r2 - 2*r1*r3 - 2*r1*r4 - 2*r2*r3 - 2*r2*r4 - 2*r3*r4 + r1*r2*r3 + r1*r2*r4 + r1*r3*r4 + r2*r3*r4 - 4));
                b4 = -(r1 + 2*r2 + 2*r3 + r4 - r2^2*r3^2 - 2*r1*r2 - 2*r1*r3 - r1*r4 - 4*r2*r3 - 2*r2*r4 - 2*r3*r4 + r1*r2^2 + r1*r3^2 + 2*r2*r3^2 + 2*r2^2*r3 + r2^2*r4 + r3^2*r4 - r2^2 - r3^2 - 2*r1*r2*r3^2 - 2*r1*r2^2*r3 - r1*r2^2*r4 - r1*r3^2*r4 - 2*r2*r3^2*r4 - 2*r2^2*r3*r4 + r1*r2^2*r3^2 + r2^2*r3^2*r4 + 4*r1*r2*r3 + 2*r1*r2*r4 + 2*r1*r3*r4 + 4*r2*r3*r4 - 4*r1*r2*r3*r4 + 2*r1*r2*r3^2*r4 + 2*r1*r2^2*r3*r4 - r1*r2^2*r3^2*r4 - 1)/((r1*r4^2 + r2*r4^2 + r3*r4^2 - r4^3 + r1*r2*r3 - r1*r2*r4 - r1*r3*r4 - r2*r3*r4)*(3*r1 + 3*r2 + 3*r3 + 3*r4 - 2*r1*r2 - 2*r1*r3 - 2*r1*r4 - 2*r2*r3 - 2*r2*r4 - 2*r3*r4 + r1*r2*r3 + r1*r2*r4 + r1*r3*r4 + r2*r3*r4 - 4));
                

                
                
            elseif fLCBDF == 1
                
                if strcmp(addTerm,'R')
                    
                    a1=(- 48*r1^4*r2*r3 - 48*r1^4*r2*r4 + 96*r1^4*r2 - 48*r1^4*r3*r4 + 96*r1^4*r3 + 96*r1^4*r4 - 144*r1^4 + 36*r1^3*r2^2*r3 + 36*r1^3*r2^2*r4 - 72*r1^3*r2^2 + 36*r1^3*r2*r3^2 + 172*r1^3*r2*r3*r4 - 124*r1^3*r2*r3 + 36*r1^3*r2*r4^2 - 124*r1^3*r2*r4 + 4*r1^3*r2 + 36*r1^3*r3^2*r4 - 72*r1^3*r3^2 + 36*r1^3*r3*r4^2 - 124*r1^3*r3*r4 + 4*r1^3*r3 - 72*r1^3*r4^2 + 4*r1^3*r4 + 188*r1^3 - 24*r1^2*r2^2*r3^2 - 99*r1^2*r2^2*r3*r4 + 39*r1^2*r2^2*r3 - 24*r1^2*r2^2*r4^2 + 39*r1^2*r2^2*r4 + 69*r1^2*r2^2 - 99*r1^2*r2*r3^2*r4 + 39*r1^2*r2*r3^2 - 99*r1^2*r2*r3*r4^2 + 78*r1^2*r2*r3*r4 + 69*r1^2*r2*r3 + 39*r1^2*r2*r4^2 + 69*r1^2*r2*r4 - 96*r1^2*r2 - 24*r1^2*r3^2*r4^2 + 39*r1^2*r3^2*r4 + 69*r1^2*r3^2 + 39*r1^2*r3*r4^2 + 69*r1^2*r3*r4 - 96*r1^2*r3 + 69*r1^2*r4^2 - 96*r1^2*r4 - 69*r1^2 + 50*r1*r2^2*r3^2*r4 - 2*r1*r2^2*r3^2 + 50*r1*r2^2*r3*r4^2 - 2*r1*r2^2*r3*r4 - 36*r1*r2^2*r3 - 2*r1*r2^2*r4^2 - 36*r1*r2^2*r4 - 22*r1*r2^2 + 50*r1*r2*r3^2*r4^2 - 2*r1*r2*r3^2*r4 - 36*r1*r2*r3^2 - 2*r1*r2*r3*r4^2 - 72*r1*r2*r3*r4 + 26*r1*r2*r3 - 36*r1*r2*r4^2 + 26*r1*r2*r4 + 46*r1*r2 - 2*r1*r3^2*r4^2 - 36*r1*r3^2*r4 - 22*r1*r3^2 - 36*r1*r3*r4^2 + 26*r1*r3*r4 + 46*r1*r3 - 22*r1*r4^2 + 46*r1*r4 - 25*r2^2*r3^2*r4^2 + r2^2*r3^2 + r2^2*r3*r4 + 11*r2^2*r3 + r2^2*r4^2 + 11*r2^2*r4 + r2*r3^2*r4 + 11*r2*r3^2 + r2*r3*r4^2 + 22*r2*r3*r4 - 23*r2*r3 + 11*r2*r4^2 - 23*r2*r4 + r3^2*r4^2 + 11*r3^2*r4 + 11*r3*r4^2 - 23*r3*r4)/(25*(r1 - r2)^2*(r1 - r3)^2*(r1 - r4)^2);
                    a2=                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  -(70*r1 + 11*r3 + 11*r4 - 34*r1*r3 - 34*r1*r4 + r3*r4 + 35*r1^2*r3 + 35*r1^2*r4 - 12*r1^3*r3 - 12*r1^3*r4 - 71*r1^2 + 24*r1^3 + r1^2*r3*r4 - 2*r1*r3*r4 - 23)/(25*(r1 - r2)*(r2 - r3)*(r1*r2 - r1*r4 + r2*r4 - r2^2));
                    a3=                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   (70*r1 + 11*r2 + 11*r4 - 34*r1*r2 - 34*r1*r4 + r2*r4 + 35*r1^2*r2 - 12*r1^3*r2 + 35*r1^2*r4 - 12*r1^3*r4 - 71*r1^2 + 24*r1^3 + r1^2*r2*r4 - 2*r1*r2*r4 - 23)/(25*(r1 - r3)*(r3 - r4)*(r1*r2 - r1*r3 - r2*r3 + r3^2));
                    a4=                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          -(70*r1 + 11*r2 + 11*r3 - 34*r1*r2 - 34*r1*r3 + r2*r3 + 35*r1^2*r2 + 35*r1^2*r3 - 12*r1^3*r2 - 12*r1^3*r3 - 71*r1^2 + 24*r1^3 + r1^2*r2*r3 - 2*r1*r2*r3 - 23)/(25*(r1 - r4)*(r1*r4^2 + r2*r4^2 + r3*r4^2 - r4^3 + r1*r2*r3 - r1*r2*r4 - r1*r3*r4 - r2*r3*r4));
                    a5=                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           -(36*r1 + 11*r2 + 11*r3 + 11*r4 - 24*r1*r2 - 24*r1*r3 - 24*r1*r4 + r2*r3 + r2*r4 + r3*r4 + 12*r1*r2*r3 + 12*r1*r2*r4 + 12*r1*r3*r4 - 13*r2*r3*r4 - 23)/(25*(r1^2*r2 + r1^2*r3 + r1^2*r4 - r1^3 - r1*r2*r3 - r1*r2*r4 - r1*r3*r4 + r2*r3*r4));
                    
                    b0 = fLC*(1-r1);
                    b1 = -(r2 - 36*r1 + r3 + r4 + 24*r1*r2 + 24*r1*r3 + 24*r1*r4 - 13*r2*r3 - 13*r2*r4 - 13*r3*r4 - 12*r1*r2*r3 - 12*r1*r2*r4 - 12*r1*r3*r4 + 25*r2*r3*r4 + 11)/(25*(r1 - r2)*(r1 - r3)*(r1 - r4));
                    b2 =                                                                                     (12*(r1 + r3 + r4 - r1*r3 - r1*r4 - r3*r4 + r1*r3*r4 - 1))/(25*(r1 - r2)*(r2 - r3)*(r2 - r4));
                    b3 =                                                                         -(12*(r1 + r2 + r4 - r1*r2 - r1*r4 - r2*r4 + r1*r2*r4 - 1))/(25*(r3 - r4)*(r1*r2 - r1*r3 - r2*r3 + r3^2));
                    b4 =                                  (12*(r1 + r2 + r3 - r1*r2 - r1*r3 - r2*r3 + r1*r2*r3 - 1))/(25*(r1*r4^2 + r2*r4^2 + r3*r4^2 - r4^3 + r1*r2*r3 - r1*r2*r4 - r1*r3*r4 - r2*r3*r4));
                    
                    
                elseif strcmp(addTerm,'L')
                    
                    a1=                                    (48*r1 + 11*r2 + 11*r3 + 11*r4 + 11*r5 - 36*r1*r2 - 36*r1*r3 - 36*r1*r4 + r2*r3 - 36*r1*r5 + r2*r4 + r2*r5 + r3*r4 + r3*r5 + r4*r5 + 24*r1*r2*r3 + 24*r1*r2*r4 + 24*r1*r2*r5 + 24*r1*r3*r4 + 24*r1*r3*r5 - 13*r2*r3*r4 + 24*r1*r4*r5 - 13*r2*r3*r5 - 13*r2*r4*r5 - 13*r3*r4*r5 - 12*r1*r2*r3*r4 - 12*r1*r2*r3*r5 - 12*r1*r2*r4*r5 - 12*r1*r3*r4*r5 + 25*r2*r3*r4*r5 - 23)/(25*(r1 - r2)*(r1 - r3)*(r1*r4 + r1*r5 - r4*r5 - r1^2));
                    a2=                                                                                                                                     -(59*r1 + 11*r3 + 11*r4 + 11*r5 - 35*r1*r3 - 35*r1*r4 - 35*r1*r5 + r3*r4 + r3*r5 + r4*r5 + 24*r1^2*r3 + 24*r1^2*r4 + 24*r1^2*r5 - 36*r1^2 - 12*r1^2*r3*r4 - 12*r1^2*r3*r5 - 12*r1^2*r4*r5 + 11*r1*r3*r4 + 11*r1*r3*r5 + 11*r1*r4*r5 - 13*r3*r4*r5 + 13*r1*r3*r4*r5 - 23)/(25*(r1 - r2)*(r2 - r3)*(r2*r4 + r2*r5 - r4*r5 - r2^2));
                    a3=                                                                                                                                     -(59*r1 + 11*r2 + 11*r4 + 11*r5 - 35*r1*r2 - 35*r1*r4 - 35*r1*r5 + r2*r4 + r2*r5 + r4*r5 + 24*r1^2*r2 + 24*r1^2*r4 + 24*r1^2*r5 - 36*r1^2 - 12*r1^2*r2*r4 - 12*r1^2*r2*r5 - 12*r1^2*r4*r5 + 11*r1*r2*r4 + 11*r1*r2*r5 + 11*r1*r4*r5 - 13*r2*r4*r5 + 13*r1*r2*r4*r5 - 23)/(25*(r3 - r4)*(r3 - r5)*(r1*r2 - r1*r3 - r2*r3 + r3^2));
                    a4=                                                                                              (59*r1 + 11*r2 + 11*r3 + 11*r5 - 35*r1*r2 - 35*r1*r3 + r2*r3 - 35*r1*r5 + r2*r5 + r3*r5 + 24*r1^2*r2 + 24*r1^2*r3 + 24*r1^2*r5 - 36*r1^2 - 12*r1^2*r2*r3 - 12*r1^2*r2*r5 - 12*r1^2*r3*r5 + 11*r1*r2*r3 + 11*r1*r2*r5 + 11*r1*r3*r5 - 13*r2*r3*r5 + 13*r1*r2*r3*r5 - 23)/(25*(r4 - r5)*(r1*r4^2 + r2*r4^2 + r3*r4^2 - r4^3 + r1*r2*r3 - r1*r2*r4 - r1*r3*r4 - r2*r3*r4));
                    a5=  (59*r1 + 11*r2 + 11*r3 + 11*r4 - 35*r1*r2 - 35*r1*r3 - 35*r1*r4 + r2*r3 + r2*r4 + r3*r4 + 24*r1^2*r2 + 24*r1^2*r3 + 24*r1^2*r4 - 36*r1^2 - 12*r1^2*r2*r3 - 12*r1^2*r2*r4 - 12*r1^2*r3*r4 + 11*r1*r2*r3 + 11*r1*r2*r4 + 11*r1*r3*r4 - 13*r2*r3*r4 + 13*r1*r2*r3*r4 - 23)/(25*(r1*r5^3 + r2*r5^3 + r3*r5^3 + r4*r5^3 - r5^4 - r1*r2*r5^2 - r1*r3*r5^2 - r1*r4*r5^2 - r2*r3*r5^2 - r2*r4*r5^2 - r3*r4*r5^2 - r1*r2*r3*r4 + r1*r2*r3*r5 + r1*r2*r4*r5 + r1*r3*r4*r5 + r2*r3*r4*r5));
                    
                    b0 = fLC*(1-r1);
                    
                    b1=                                         -(12*(r2 + r3 + r4 - r2*r3 - r2*r4 - r3*r4 + r2*r3*r4 - 1))/(25*(r1 - r2)*(r1 - r3)*(r1 - r4));
                    b2=                                          (12*(r1 + r3 + r4 - r1*r3 - r1*r4 - r3*r4 + r1*r3*r4 - 1))/(25*(r1 - r2)*(r2 - r3)*(r2 - r4));
                    b3=                              -(12*(r1 + r2 + r4 - r1*r2 - r1*r4 - r2*r4 + r1*r2*r4 - 1))/(25*(r3 - r4)*(r1*r2 - r1*r3 - r2*r3 + r3^2));
                    b4=(12*(r1 + r2 + r3 - r1*r2 - r1*r3 - r2*r3 + r1*r2*r3 - 1))/(25*(r1*r4^2 + r2*r4^2 + r3*r4^2 - r4^3 + r1*r2*r3 - r1*r2*r4 - r1*r3*r4 - r2*r3*r4));
                    
                end
            end
        end
    end
    
    for m = 1:length(xv)
        
        for n = 1:length(yv)
            
            y = 1i*yv(n);
            x =    xv(m);
            
            if order == 4
                
                if fixStep == 1
                    
                    exVE = [25/12-x    -(4*y+4) (6*y+3) -(4*y+4/3) 1/4+y];
                    
                else
                    
                    x = x/(1-r1);
                    
                    if fLCBDF == 0
                        
                        exVE = [1-b0*x a1-b1*y a2-b2*y a3-b3*y a4-b4*y];
                        
                    elseif fLCBDF == 1
                        
                        if strcmp(addTerm,'R')
                            
                            exVE = [1-b0*x a1-b1*y-a5*x a2-b2*y a3-b3*y a4-b4*y];
                            
                        elseif strcmp(addTerm,'L')
                            
                            exVE = [1-b0*x a1-b1*y a2-b2*y a3-b3*y a4-b4*y a5];
                            
                        end
                        
                    end
                    
                    
                end
                
            elseif order == 2
                
                if fixStep == 1
                    
                    exVE = [1-2/3*x -4/3-4/3*y 1/3+2/3*y];
                    
                else
                    
                    x = x/(1-r1);
                    exVE = [1-b0*x a1-b1*y a2-b2*y];
                    
                end
                
            end
            
            
            pRoots  = roots(exVE);
            p1(n,m) = max(abs(pRoots));
            
            
        end
    end
    
    %h = figure(1);
    contour(X,Y,p1,[1 1],'LineWidth',2,'LineColor',lC(k)); %contour for |R(z)|=1
    grid on
    %,'LineColor','r'
    hold on
    
end

% legendStr = [];
% for i = 1:5
%     r = 1 + (i-3)*0.2;
%     legendStr = [legendStr, num2str(r)];
% end
%
%
% legend(legendStr);
titleStr = ['newTimeLevel',num2str(newTimeLevel)];
title(titleStr);
%print(titleStr,'-dpng')


end