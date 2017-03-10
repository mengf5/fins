function [L1] = formLHSVelU1New(fS)

Nxg = fS.Nxg;
Nyg = fS.Nyg;

Nx = fS.Nx;
Ny = fS.Ny;

gL = fS.ghostLine;

hx  = fS.hx;
hy  = fS.hy;

BC  = fS.BC; % BC ==4 means dirichlet boundary conditions with exact ghost points 
             % corner points are still calculated with Taylor expansion 

             % BC ==6 periodic in x direction, Dirichlet in y direction 
             
extOrder   = fS.extOrder;       
imTime     = fS.imTime;
mu         = fS.mu; 

%% get interior pts index 
M    = Nyg*Nxg;

L1 = spalloc(M,M,7*M);

[intPts,lIntPts] = getInteriorIndex(fS);

coeffC = ones(1,lIntPts);%*(1-imTime(1)*mu*(-(30/(12*hy^2)+30/(12*hx^2))));
L1 = L1 + sparse(intPts,intPts,coeffC,M,M);
%
L1 = L1 + sparse(intPts,intPts+1,coeffC,M,M);
L1 = L1 + sparse(intPts,intPts-1,coeffC,M,M);
L1 = L1 + sparse(intPts,intPts+2,coeffC,M,M);
L1 = L1 + sparse(intPts,intPts-2,coeffC,M,M);
%
L1 = L1 + sparse(intPts,intPts+Nyg,coeffC,M,M);
L1 = L1 + sparse(intPts,intPts-Nyg,coeffC,M,M);
L1 = L1 + sparse(intPts,intPts+2*Nyg,coeffC,M,M);
L1 = L1 + sparse(intPts,intPts-2*Nyg,coeffC,M,M);
%% 
%         L1(ix,ix +1) =   -imTime(1)*mu*16/(12*hy^2);
%         L1(ix,ix -1) =   -imTime(1)*mu*16/(12*hy^2);
%         L1(ix,ix +2) =   -imTime(1)*mu*(-1/(12*hy^2));
%         L1(ix,ix -2) =   -imTime(1)*mu*(-1/(12*hy^2));
%         
%         L1(ix,ix + Nyg)   =   -imTime(1)*mu*16/(12*hx^2);
%         L1(ix,ix - Nyg)   =   -imTime(1)*mu*16/(12*hx^2);
%         L1(ix,ix + 2*Nyg) =   -imTime(1)*mu*(-1/(12*hx^2));
%         L1(ix,ix - 2*Nyg) =   -imTime(1)*mu*(-1/(12*hx^2));

axis = 0;
side = 0;
pos =  1;
%localBC = BC(axis,side);


%% evaluate the boundary line, ghost line indexing for given (axis,side,position)
%  position 0: boundary 1:first ghost line 2:second ghost line

bcxStart = (side==0)*ia + (side==1)*ib;
bcyStart = (side==0)*ja + (side==1)*jb;

bcxEnd = (axis==1)*((side==1)*ia + (side==0)*ib) +  (axis==0)*(bcxStart);
bcyEnd = (axis==0)*((side==1)*ja + (side==0)*jb) +  (axis==1)*(bcyStart);

lBcx = abs(bcxEnd - bcxStart) + 1;
lBcy = abs(bcyEnd - bcyStart) + 1;

bcx = bcxStart:bcxEnd;
bcy = bcyStart:bcyEnd;

bcx = (axis==0)*(bcx - (-1)^side * pos );
bcy = (axis==1)*(bcy - (-1)^side * pos );

bcPts = kron((bcx-1)*Nyg,ones(1,lBcx)) + kron(ones(1,lBcy),bcy);
lBcPts= length(bcPts);
%% 

coeffB = ones(1,length(lBcPts));
L1     = L1 + sparse(bcPts,bcPts,coeffB,M,M);

    






for ix = 1:3
    
    if BC == 2 || BC == 6
        
        L1(ix,ix        )     =   1;
        L1(ix,ix + M - 5*Nyg) =  -1;
        
    else
        
        L1(ix,ix) = 1;
    end
    
end

for ix = 4:Nyg-3

    if BC == 1
                
%         L1(ix,ix        ) =  -1/(12*hx^2);
%         L1(ix,ix +   Nyg) =  16/(12*hx^2);
%         L1(ix,ix + 2*Nyg) = -30/(12*hx^2);
%         L1(ix,ix + 3*Nyg) =  16/(12*hx^2);
%         L1(ix,ix + 4*Nyg) =  -1/(12*hx^2);
        
        L1(ix,ix        ) =  -1/(2*hx^3);
        L1(ix,ix +   Nyg) =   1/(  hx^3);
        L1(ix,ix + 3*Nyg) =  -1/(  hx^3);
        L1(ix,ix + 4*Nyg) =   1/(2*hx^3);
        
    
    elseif BC == 2 || BC == 6
        
        L1(ix,ix        )     =   1;
        L1(ix,ix + M - 5*Nyg) =  -1;
        
    elseif BC == 4
        L1(ix,ix) = 1;
    end
    
end

for ix = Nyg-2:Nyg+3
    
    if BC == 2 || BC == 6
        
        L1(ix,ix        )     =   1;
        L1(ix,ix + M - 5*Nyg) =  -1;
        
    else
        
        L1(ix,ix) = 1;
    end
    
end

if BC == 1
    %% fixed corner
    ix = Nyg+2;
    L1(ix,ix) = 1;
    L1(ix,ix +   Nyg + 1) = -15/4;
    L1(ix,ix + 2*Nyg + 2) = 3;
    L1(ix,ix + 3*Nyg + 3) = -1/4;
    
    ix = Nyg+1;
    L1(ix,ix) = 1;
    L1(ix,ix +   Nyg + 1*2) = -15/4;
    L1(ix,ix + 2*Nyg + 2*2) = 3;
    L1(ix,ix + 3*Nyg + 3*2) = -1/4;
    
    ix = 2;
    L1(ix,ix) = 1;
    L1(ix,ix + 2*  Nyg + 1) = -15/4;
    L1(ix,ix + 2*2*Nyg + 2) = 3;
    L1(ix,ix + 2*3*Nyg + 3) = -1/4;
end

for ix = Nyg + 4: Nyg + Nyg-3
    
    if BC == 1
        
        L1(ix,ix -   Nyg) =  ( 1)/(12*hx);
        L1(ix,ix        ) =  (-8)/(12*hx);
        L1(ix,ix + 2*Nyg) =  ( 8)/(12*hx);
        L1(ix,ix + 3*Nyg) =  (-1)/(12*hx);
    
    elseif BC == 2 || BC == 6
        
        L1(ix,ix        )     =   1;
        L1(ix,ix + M - 5*Nyg) =  -1;
        
    elseif BC == 4
        L1(ix,ix) = 1;
    end
    
end

for ix = 2*Nyg-2:2*Nyg
    
    if BC == 2 || BC == 6
        
        L1(ix,ix        )     =   1;
        L1(ix,ix + M - 5*Nyg) =  -1;
        
    else
        
        L1(ix,ix) = 1;
    end
    
end

% fix corner 
if BC == 1
    ix = 2*Nyg-1;
    L1(ix,ix) = 1;
    L1(ix,ix +   Nyg - 1) = -15/4;
    L1(ix,ix + 2*Nyg - 2) = 3;
    L1(ix,ix + 3*Nyg - 3) = -1/4;
end

if BC == 2 || BC == 6
    
    indexInterior = 2*Nyg+1;
    
else
    
    for ix = 2*Nyg+1:3*Nyg        
        L1(ix,ix) = 1;
    end
    
    indexInterior = 3*Nyg+1;
    
end

for ix = indexInterior:M-3*Nyg
    
    if mod(ix - 2*Nyg,Nyg) == 1
        
        if BC == 1 || BC == 6
            
            if extOrder==6
                
                L1(ix,ix    ) =   1;
                L1(ix,ix + 1) =  -6;
                L1(ix,ix + 2) =  15;
                L1(ix,ix + 3) = -20;
                L1(ix,ix + 4) =  15;
                L1(ix,ix + 5) =  -6;
                L1(ix,ix + 6) =   1;
                
            elseif extOrder==5
                
                L1(ix,ix    ) =   1;
                L1(ix,ix + 1) =  -5;
                L1(ix,ix + 2) =  10;
                L1(ix,ix + 3) = -10;
                L1(ix,ix + 4) =   5;
                L1(ix,ix + 5) =  -1;
                
            elseif extOrder==4
                
                L1(ix,ix    ) =   1;
                L1(ix,ix + 1) =  -4;
                L1(ix,ix + 2) =   6;
                L1(ix,ix + 3) =  -4;
                L1(ix,ix + 4) =   1;
                
            end
            
        elseif BC == 2
            
            L1(ix,ix)            = 1;
            L1(ix,ix + Nyg - 5 ) = -1;
            
        elseif BC == 4
            L1(ix,ix) = 1;
        end
        
    elseif mod(ix - 2*Nyg,Nyg) == 0
        
        if BC == 1 || BC == 6
            
            if extOrder==6
                
                L1(ix,ix    ) =   1;
                L1(ix,ix - 1) =  -6;
                L1(ix,ix - 2) =  15;
                L1(ix,ix - 3) = -20;
                L1(ix,ix - 4) =  15;
                L1(ix,ix - 5) =  -6;
                L1(ix,ix - 6) =   1;
            elseif extOrder==5
                
                L1(ix,ix    ) =   1;
                L1(ix,ix - 1) =  -5;
                L1(ix,ix - 2) =  10;
                L1(ix,ix - 3) = -10;
                L1(ix,ix - 4) =   5;
                L1(ix,ix - 5) =  -1;
                
            elseif extOrder==4
                
                L1(ix,ix    ) =   1;
                L1(ix,ix - 1) =  -4;
                L1(ix,ix - 2) =   6;
                L1(ix,ix - 3) =  -4;
                L1(ix,ix - 4) =   1;                
                
            end
            
        elseif BC == 2
            L1(ix,ix)            = 1;
            L1(ix,ix - Nyg + 5 ) = -1;
            
        elseif BC == 4
            L1(ix,ix) = 1;
        end
        
    elseif mod(ix - 2*Nyg,Nyg) ==  2
        
        if BC == 1 || BC == 6
            
            L1(ix,ix - 1) =   -1/(12*hy^2);
            L1(ix,ix    ) =   16/(12*hy^2);
            L1(ix,ix + 1) =  -30/(12*hy^2);
            L1(ix,ix + 2) =   16/(12*hy^2);
            L1(ix,ix + 3) =   -1/(12*hy^2);
            
        elseif BC == 2
            L1(ix,ix)            = 1;
            L1(ix,ix + Nyg - 5 ) = -1;
            
        elseif BC == 4
            L1(ix,ix) = 1;
        end
        
        
    elseif mod(ix - 2*Nyg,Nyg) ==  Nyg-1

        if BC == 1 || BC == 6
            
            L1(ix,ix + 1) =   -1/(12*hy^2);
            L1(ix,ix    ) =   16/(12*hy^2);
            L1(ix,ix - 1) =  -30/(12*hy^2);
            L1(ix,ix - 2) =   16/(12*hy^2);
            L1(ix,ix - 3) =   -1/(12*hy^2);
        
        elseif BC == 2
            L1(ix,ix)            = 1;
            L1(ix,ix - Nyg + 5 ) = -1;
            
        elseif BC == 4
            L1(ix,ix) = 1;
        end
        
    elseif mod(ix - 2*Nyg,Nyg) ==  3
        
        if BC == 2
            L1(ix,ix   ) =  1-imTime(1)*mu*(-(30/(12*hy^2)+30/(12*hx^2)));
            L1(ix,ix +1) =   -imTime(1)*mu*16/(12*hy^2);
            L1(ix,ix -1) =   -imTime(1)*mu*16/(12*hy^2);
            L1(ix,ix +2) =   -imTime(1)*mu*(-1/(12*hy^2));
            L1(ix,ix -2) =   -imTime(1)*mu*(-1/(12*hy^2));
            
            L1(ix,ix + Nyg)   =   -imTime(1)*mu*16/(12*hx^2);
            L1(ix,ix - Nyg)   =   -imTime(1)*mu*16/(12*hx^2);
            L1(ix,ix + 2*Nyg) =   -imTime(1)*mu*(-1/(12*hx^2));
            L1(ix,ix - 2*Nyg) =   -imTime(1)*mu*(-1/(12*hx^2));
        else
            
            L1(ix,ix)     =   1;
            
        end
        
    elseif mod(ix - 2*Nyg,Nyg) ==  Nyg-2
        
        
        if BC == 2
            
            L1(ix,ix)               =   1;
            L1(ix,ix - Nyg + 5)     =   -1;
            
        else
            
            L1(ix,ix)     =   1;
            
        end
        
    else
        L1(ix,ix   ) =  1-imTime(1)*mu*(-(30/(12*hy^2)+30/(12*hx^2)));
        L1(ix,ix +1) =   -imTime(1)*mu*16/(12*hy^2);
        L1(ix,ix -1) =   -imTime(1)*mu*16/(12*hy^2);
        L1(ix,ix +2) =   -imTime(1)*mu*(-1/(12*hy^2));
        L1(ix,ix -2) =   -imTime(1)*mu*(-1/(12*hy^2));
        
        L1(ix,ix + Nyg)   =   -imTime(1)*mu*16/(12*hx^2);
        L1(ix,ix - Nyg)   =   -imTime(1)*mu*16/(12*hx^2);
        L1(ix,ix + 2*Nyg) =   -imTime(1)*mu*(-1/(12*hx^2));
        L1(ix,ix - 2*Nyg) =   -imTime(1)*mu*(-1/(12*hx^2));
        
    end
    
end

for ix = M-3*Nyg + 1 : M-2*Nyg
    
    if BC == 2 || BC == 6
    
        L1(ix,ix) = 1;
        L1(ix,ix - M + 5*Nyg) = -1;
        
    else
        
    L1(ix,ix) = 1;
    
    end
end

for ix = M-2*Nyg + 1 : M-2*Nyg + 3

    if BC == 2 || BC == 6
    
        L1(ix,ix) = 1;
        L1(ix,ix - M + 5*Nyg) = -1;
        
    else
        
    L1(ix,ix) = 1;
    
    end
    
end

for ix = M-2*Nyg + 4 : M-Nyg-3
    
    if BC == 1
        
        L1(ix,ix        ) =  ( 8)/(12*hx);
        L1(ix,ix +   Nyg) =  (-1)/(12*hx);
        L1(ix,ix - 2*Nyg) =  (-8)/(12*hx);
        L1(ix,ix - 3*Nyg) =  ( 1)/(12*hx);
    
    elseif BC == 2 || BC == 6
        
        L1(ix,ix) = 1;
        L1(ix,ix - M + 5*Nyg) = -1;
        
    elseif BC == 4
        
        L1(ix,ix) = 1;
    
    end
    
end

for ix = M-Nyg-2 : M - Nyg + 3
    
    if BC == 2 || BC == 6
    
        L1(ix,ix) = 1;
        L1(ix,ix - M + 5*Nyg) = -1;
        
    else
        
    L1(ix,ix) = 1;
    
    end
end


for ix = M - Nyg + 4: M-3
    
    if BC == 1
        
        L1(ix,ix        ) =  -1/(12*hx^2);
        L1(ix,ix -   Nyg) =  16/(12*hx^2);
        L1(ix,ix - 2*Nyg) = -30/(12*hx^2);
        L1(ix,ix - 3*Nyg) =  16/(12*hx^2);
        L1(ix,ix - 4*Nyg) =  -1/(12*hx^2);
                
    elseif BC == 2 || BC == 6
        
        L1(ix,ix) = 1;
        L1(ix,ix - M + 5*Nyg) = -1;
        
    elseif BC == 4
        L1(ix,ix) = 1;
    end
end

for ix = M-2:M
    
    if BC == 2 || BC == 6
    
        L1(ix,ix) = 1;
        L1(ix,ix - M + 5*Nyg) = -1;
        
    else
        
    L1(ix,ix) = 1;
    
    end
    
end

L1 = sparse(L1);
end


function setupBCs(dir,L1)




end
