function [L1] = formLHSVelV1(fS)

Nxg = fS.Nxg;
Nyg = fS.Nyg;

hx  = fS.hx;
hy  = fS.hy;

BC  = fS.BC; % BC ==4 means dirichlet boundary conditions with exact ghost points
% corner points are still calculated with Taylor expansion

extOrder   = fS.extOrder;
imTime  = fS.imTime;
mu = fS.mu;

M    = Nyg*Nxg;

L1 = spalloc(M,M,7*M);
%L1   = zeros(M,M);

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
        
        if extOrder==6
            L1(ix,ix        ) =   1;
            L1(ix,ix +   Nyg) =  -6;
            L1(ix,ix + 2*Nyg) =  15;
            L1(ix,ix + 3*Nyg) = -20;
            L1(ix,ix + 4*Nyg) =  15;
            L1(ix,ix + 5*Nyg) =  -6;
            L1(ix,ix + 6*Nyg) =   1;
        elseif extOrder==5
            L1(ix,ix        ) =   1;
            L1(ix,ix +   Nyg) =  -5;
            L1(ix,ix + 2*Nyg) =  10;
            L1(ix,ix + 3*Nyg) = -10;
            L1(ix,ix + 4*Nyg) =   5;
            L1(ix,ix + 5*Nyg) =  -1;
            
        elseif extOrder==4
            L1(ix,ix        ) =   1;
            L1(ix,ix +   Nyg) =  -4;
            L1(ix,ix + 2*Nyg) =   6;
            L1(ix,ix + 3*Nyg) =  -4;
            L1(ix,ix + 4*Nyg) =   1;
            
        end
    
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
% fix corner 
if BC == 1
    
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
        
        L1(ix,ix -   Nyg) =  ( -1)/(12*hx^2);
        L1(ix,ix        ) =  ( 16)/(12*hx^2);
        L1(ix,ix +   Nyg) =  (-30)/(12*hx^2);
        L1(ix,ix + 2*Nyg) =  ( 16)/(12*hx^2);
        L1(ix,ix + 3*Nyg) =  ( -1)/(12*hx^2);
        
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
            
            L1(ix,ix    ) =   -1/(12*hy^2);
            L1(ix,ix + 1) =   16/(12*hy^2);
            L1(ix,ix + 2) =  -30/(12*hy^2);
            L1(ix,ix + 3) =   16/(12*hy^2);
            L1(ix,ix + 4) =   -1/(12*hy^2);

        elseif BC == 2
            L1(ix,ix)            = 1;
            L1(ix,ix + Nyg - 5 ) = -1;
            
        elseif BC == 4
            L1(ix,ix) = 1;
        end
        
    elseif mod(ix - 2*Nyg,Nyg) == 0
        
        if BC == 1 || BC == 6
            
            L1(ix,ix    ) =   -1/(12*hy^2);
            L1(ix,ix - 1) =   16/(12*hy^2);
            L1(ix,ix - 2) =  -30/(12*hy^2);
            L1(ix,ix - 3) =   16/(12*hy^2);
            L1(ix,ix - 4) =   -1/(12*hy^2);
        
        elseif BC == 2
            
            L1(ix,ix)            = 1;
            L1(ix,ix - Nyg + 5 ) = -1;
            
        elseif BC == 4
            L1(ix,ix) = 1;
        end
        
    elseif mod(ix - 2*Nyg,Nyg) ==  2
        
        if BC == 1 || BC == 6
            
            L1(ix,ix - 1) =   1/(12*hy);
            L1(ix,ix    ) =  -8/(12*hy);
            L1(ix,ix + 2) =   8/(12*hy);
            L1(ix,ix + 3) =  -1/(12*hy);
            
        elseif BC == 2
            L1(ix,ix)            = 1;
            L1(ix,ix + Nyg - 5 ) = -1;
            
        elseif BC == 4
            L1(ix,ix) = 1;
        end
        
    elseif mod(ix - 2*Nyg,Nyg) ==  Nyg-1
        
        if BC == 1 || BC == 6
            
            L1(ix,ix + 1) =  -1/(12*hy);
            L1(ix,ix    ) =   8/(12*hy);
            L1(ix,ix - 2) =  -8/(12*hy);
            L1(ix,ix - 3) =   1/(12*hy);
        
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
        
        L1(ix,ix +   Nyg) =  ( -1)/(12*hx^2);
        L1(ix,ix        ) =  ( 16)/(12*hx^2);
        L1(ix,ix -   Nyg) =  (-30)/(12*hx^2);
        L1(ix,ix - 2*Nyg) =  ( 16)/(12*hx^2);
        L1(ix,ix - 3*Nyg) =  ( -1)/(12*hx^2);
        
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
        
        if extOrder==6
            L1(ix,ix        ) =   1;
            L1(ix,ix -   Nyg) =  -6;
            L1(ix,ix - 2*Nyg) =  15;
            L1(ix,ix - 3*Nyg) = -20;
            L1(ix,ix - 4*Nyg) =  15;
            L1(ix,ix - 5*Nyg) =  -6;
            L1(ix,ix - 6*Nyg) =   1;
            
        elseif extOrder==5
            L1(ix,ix        ) =   1;
            L1(ix,ix -   Nyg) =  -5;
            L1(ix,ix - 2*Nyg) =  10;
            L1(ix,ix - 3*Nyg) = -10;
            L1(ix,ix - 4*Nyg) =   5;
            L1(ix,ix - 5*Nyg) =  -1;
            
        elseif extOrder==4
            L1(ix,ix        ) =   1;
            L1(ix,ix -   Nyg) =  -4;
            L1(ix,ix - 2*Nyg) =   6;
            L1(ix,ix - 3*Nyg) =  -4;
            L1(ix,ix - 4*Nyg) =   1;
            
        end

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






