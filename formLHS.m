function [lhsU] = formLHS(fS)

% get info from INS class 
Nxg = fS.Nxg;
Nyg = fS.Nyg;

hx  = fS.hx;
hy  = fS.hy;

BC  = fS.BC;  
             
extOrder   = fS.extOrder;       
imTime     = fS.imTime;
mu         = fS.mu; 

% total number of grid points(including ghosts)
M    = Nyg*Nxg;
% allocate LHS matrix for U
lhsU = spalloc(M,M,7*M);
lhsV = spalloc(M,M,7*M);

% get the indecies for all the interior points for U
[intPts,lIntPts] = getIndex(fS,ia+1,ib-1,ja+1,jb-1);
% put the laplacian operator onto these indecies
lhsU = setLaplacian(lhsU,hx,hy,mu,imTime,M,intPts,lIntPts);
lhsV = setLaplacian(lhsV,hx,hy,mu,imTime,M,intPts,lIntPts);

% assembling boundary conditions
% x-direction 
for axis = 0:1
for side = 0:1

  localBC = BC(axis,side);
  % on the boundary
  pos = 0;
  [bcPts,lBcPts] = getBCGLIndex(fS,axis,side,pos);
  
  switch localBC
	 
    case 1
      % u = u(t)
      pm    = 1;
      lhsU  = setOne(lhsU,M,pm,bcPts,lBcPts);
      lhsV  = setOne(lhsV,M,pm,bcPts,lBcPts);

    case 2
      % Laplace(U) = ...
      lhsU  = setLaplacian(lhsU,hx,hy,mu,imTime,M,bcPts,lBcPts);
      lhsV  = setLaplacian(lhsV,hx,hy,mu,imTime,M,bcPts,lBcPts);

    case 3

      % u = u(t)
      pm    = 1;
      lhsU  = setOne(lhsU,M,pm,bcPts,lBcPts);
      lhsV  = setOne(lhsV,M,pm,bcPts,lBcPts);
      
    case 4

    case 5
      
  end

  % on the ghost lines 

  switch localBC
	 
    case 1

      % u_x = -v_y

      pos = 1;
      
      [bcPts,lBcPts] = getBCGLIndex(fS,axis,side,pos);

      if axis == 0 
	lhsU  = setDx(lhsU,hx,side,Nyg,M,bcPts,lBcPts)
	lhsV  = setDxx(lhsV,hx,side,Nyg,M,bcPts,lBcPts)

      elseif axis == 1
	lhsU  = setDxx(lhsU,hy,side,1,M,bcPts,lBcPts)
	lhsV  = setDx(lhsU,hy,side,1,M,bcPts,lBcPts)

      end

      pos = 2;
      
      [bcPts,lBcPts] = getBCGLIndex(fS,axis,side,pos);

      if axis == 0 
	lhsU  = setDxxx(lhsU,hx,side,Nyg,M,bcPts,lBcPts)
	lhsV  = setExt(lhsV,hx,side,Nyg,M,bcPts,lBcPts)

      elseif axis == 1
	lhsU  = setExt(lhsU,hy,side,1,M,bcPts,lBcPts)
	lhsV  = setDxx(lhsU,hy,side,1,M,bcPts,lBcPts)
		
      end




    case 2
      % u(ia-i,j)  = u(Nx - i,j)

      for pos = 1:2;
      
      [bcPts,lBcPts] = getBCGLIndex(fS,axis,side,pos);

      lhsU  = setOne(lhsU,M,+1,bcPts,lBcPts);
      lhsV  = setOne(lhsV,M,+1,bcPts,lBcPts);
      
      [matchPts,lMatchPts] =  getBCGLIndex(fS,axis,side,-pos);
      
      lhsU  = setOne(lhsU,M,-1,matchPts,lMatchPts);
      lhsV  = setOne(lhsV,M,-1,matchPts,lMatchPts);

      end
      
    case 3
      % u = u(t)
      for pos = 1:2;
      
      [bcPts,lBcPts] = getBCGLIndex(fS,axis,side,pos);

      pm    = 1;
      lhsU  = setOne(lhsU,M,pm,bcPts,lBcPts);
      lhsV  = setOne(lhsV,M,pm,bcPts,lBcPts);

      end
      
    case 4

    case 5
      
  end
  
end  
end


% fix the ghost corners


    for sideX = 0:1
      for sideY = 0:1

	px = (sideX==0)*ia + (sideX==1)*ib;
	py = (sideY==0)*ja + (sideY==1)*jb;

	px1 = px - (-1)^sideX * 1;
	py1 = py - (-1)^sideY * 1;
	
	px2 = px - (-1)^sideX * 2;
	py2 = py - (-1)^sideY * 1;

	px3 = px - (-1)^sideX * 1;
	py3 = py - (-1)^sideY * 2;

	px4 = px - (-1)^sideX * 2;
	py4 = py - (-1)^sideY * 2;

	p1Index = getIndex(fS,px1,px1,py1,py1);
	p2Index = getIndex(fS,px1,px1,py1,py1);
	p3Index = getIndex(fS,px1,px1,py1,py1);
	p4Index = getIndex(fS,px1,px1,py1,py1);


	p1ExtIndex = getIndex(fS,px1,px1,py1,py1);
	
	switch localBC
	       
	  case 1
	    
	    px11 = px1 + (-1)^(sidex)*1
	    py11 = 


	    lhsU(ix,ix +   Nyg + 1) = -15/4;
	    lhsU(ix,ix + 2*Nyg + 2) = 3;
	    lhsU(ix,ix + 3*Nyg + 3) = -1/4;

	    coeff = pm*ones(1,length(lPts));
	    L     = L + sparse(pts,pts,coeff,M,M);


	  case 2
	    
	  case 3


	end

	
      end
    end    
    
    

    


if BC == 1
    %% fixed corner
    ix = Nyg+2;
    lhsU(ix,ix) = 1;
    lhsU(ix,ix +   Nyg + 1) = -15/4;
    lhsU(ix,ix + 2*Nyg + 2) = 3;
    lhsU(ix,ix + 3*Nyg + 3) = -1/4;
    
    ix = Nyg+1;
    lhsU(ix,ix) = 1;
    lhsU(ix,ix +   Nyg + 1*2) = -15/4;
    lhsU(ix,ix + 2*Nyg + 2*2) = 3;
    lhsU(ix,ix + 3*Nyg + 3*2) = -1/4;
    
    ix = 2;
    lhsU(ix,ix) = 1;
    lhsU(ix,ix + 2*  Nyg + 1) = -15/4;
    lhsU(ix,ix + 2*2*Nyg + 2) = 3;
    lhsU(ix,ix + 2*3*Nyg + 3) = -1/4;
end

for ix = Nyg + 4: Nyg + Nyg-3
    
    if BC == 1
        
        lhsU(ix,ix -   Nyg) =  ( 1)/(12*hx);
        lhsU(ix,ix        ) =  (-8)/(12*hx);
        lhsU(ix,ix + 2*Nyg) =  ( 8)/(12*hx);
        lhsU(ix,ix + 3*Nyg) =  (-1)/(12*hx);
    
    elseif BC == 2 || BC == 6
        
        lhsU(ix,ix        )     =   1;
        lhsU(ix,ix + M - 5*Nyg) =  -1;
        
    elseif BC == 4
        lhsU(ix,ix) = 1;
    end
    
end

for ix = 2*Nyg-2:2*Nyg
    
    if BC == 2 || BC == 6
        
        lhsU(ix,ix        )     =   1;
        lhsU(ix,ix + M - 5*Nyg) =  -1;
        
    else
        
        lhsU(ix,ix) = 1;
    end
    
end

% fix corner 
if BC == 1
    ix = 2*Nyg-1;
    lhsU(ix,ix) = 1;
    lhsU(ix,ix +   Nyg - 1) = -15/4;
    lhsU(ix,ix + 2*Nyg - 2) = 3;
    lhsU(ix,ix + 3*Nyg - 3) = -1/4;
end

if BC == 2 || BC == 6
    
    indexInterior = 2*Nyg+1;
    
else
    
    for ix = 2*Nyg+1:3*Nyg        
        lhsU(ix,ix) = 1;
    end
    
    indexInterior = 3*Nyg+1;
    
end

for ix = indexInterior:M-3*Nyg
    
    if mod(ix - 2*Nyg,Nyg) == 1
        
        if BC == 1 || BC == 6
            
            if extOrder==6
                
                lhsU(ix,ix    ) =   1;
                lhsU(ix,ix + 1) =  -6;
                lhsU(ix,ix + 2) =  15;
                lhsU(ix,ix + 3) = -20;
                lhsU(ix,ix + 4) =  15;
                lhsU(ix,ix + 5) =  -6;
                lhsU(ix,ix + 6) =   1;
                
            elseif extOrder==5
                
                lhsU(ix,ix    ) =   1;
                lhsU(ix,ix + 1) =  -5;
                lhsU(ix,ix + 2) =  10;
                lhsU(ix,ix + 3) = -10;
                lhsU(ix,ix + 4) =   5;
                lhsU(ix,ix + 5) =  -1;
                
            elseif extOrder==4
                
                lhsU(ix,ix    ) =   1;
                lhsU(ix,ix + 1) =  -4;
                lhsU(ix,ix + 2) =   6;
                lhsU(ix,ix + 3) =  -4;
                lhsU(ix,ix + 4) =   1;
                
            end
            
        elseif BC == 2
            
            lhsU(ix,ix)            = 1;
            lhsU(ix,ix + Nyg - 5 ) = -1;
            
        elseif BC == 4
            lhsU(ix,ix) = 1;
        end
        
    elseif mod(ix - 2*Nyg,Nyg) == 0
        
        if BC == 1 || BC == 6
            
            if extOrder==6
                
                lhsU(ix,ix    ) =   1;
                lhsU(ix,ix - 1) =  -6;
                lhsU(ix,ix - 2) =  15;
                lhsU(ix,ix - 3) = -20;
                lhsU(ix,ix - 4) =  15;
                lhsU(ix,ix - 5) =  -6;
                lhsU(ix,ix - 6) =   1;
            elseif extOrder==5
                
                lhsU(ix,ix    ) =   1;
                lhsU(ix,ix - 1) =  -5;
                lhsU(ix,ix - 2) =  10;
                lhsU(ix,ix - 3) = -10;
                lhsU(ix,ix - 4) =   5;
                lhsU(ix,ix - 5) =  -1;
                
            elseif extOrder==4
                
                lhsU(ix,ix    ) =   1;
                lhsU(ix,ix - 1) =  -4;
                lhsU(ix,ix - 2) =   6;
                lhsU(ix,ix - 3) =  -4;
                lhsU(ix,ix - 4) =   1;                
                
            end
            
        elseif BC == 2
            lhsU(ix,ix)            = 1;
            lhsU(ix,ix - Nyg + 5 ) = -1;
            
        elseif BC == 4
            lhsU(ix,ix) = 1;
        end
        
    elseif mod(ix - 2*Nyg,Nyg) ==  2
        
        if BC == 1 || BC == 6
            
            lhsU(ix,ix - 1) =   -1/(12*hy^2);
            lhsU(ix,ix    ) =   16/(12*hy^2);
            lhsU(ix,ix + 1) =  -30/(12*hy^2);
            lhsU(ix,ix + 2) =   16/(12*hy^2);
            lhsU(ix,ix + 3) =   -1/(12*hy^2);
            
        elseif BC == 2
            lhsU(ix,ix)            = 1;
            lhsU(ix,ix + Nyg - 5 ) = -1;
            
        elseif BC == 4
            lhsU(ix,ix) = 1;
        end
        
        
    elseif mod(ix - 2*Nyg,Nyg) ==  Nyg-1

        if BC == 1 || BC == 6
            
            lhsU(ix,ix + 1) =   -1/(12*hy^2);
            lhsU(ix,ix    ) =   16/(12*hy^2);
            lhsU(ix,ix - 1) =  -30/(12*hy^2);
            lhsU(ix,ix - 2) =   16/(12*hy^2);
            lhsU(ix,ix - 3) =   -1/(12*hy^2);
        
        elseif BC == 2
            lhsU(ix,ix)            = 1;
            lhsU(ix,ix - Nyg + 5 ) = -1;
            
        elseif BC == 4
            lhsU(ix,ix) = 1;
        end
        
    elseif mod(ix - 2*Nyg,Nyg) ==  3
        
        if BC == 2
            lhsU(ix,ix   ) =  1-imTime(1)*mu*(-(30/(12*hy^2)+30/(12*hx^2)));
            lhsU(ix,ix +1) =   -imTime(1)*mu*16/(12*hy^2);
            lhsU(ix,ix -1) =   -imTime(1)*mu*16/(12*hy^2);
            lhsU(ix,ix +2) =   -imTime(1)*mu*(-1/(12*hy^2));
            lhsU(ix,ix -2) =   -imTime(1)*mu*(-1/(12*hy^2));
            
            lhsU(ix,ix + Nyg)   =   -imTime(1)*mu*16/(12*hx^2);
            lhsU(ix,ix - Nyg)   =   -imTime(1)*mu*16/(12*hx^2);
            lhsU(ix,ix + 2*Nyg) =   -imTime(1)*mu*(-1/(12*hx^2));
            lhsU(ix,ix - 2*Nyg) =   -imTime(1)*mu*(-1/(12*hx^2));
        else
            
            lhsU(ix,ix)     =   1;
            
        end
        
    elseif mod(ix - 2*Nyg,Nyg) ==  Nyg-2
        
        
        if BC == 2
            
            lhsU(ix,ix)               =   1;
            lhsU(ix,ix - Nyg + 5)     =   -1;
            
        else
            
            lhsU(ix,ix)     =   1;
            
        end
        
    else
        lhsU(ix,ix   ) =  1-imTime(1)*mu*(-(30/(12*hy^2)+30/(12*hx^2)));
        lhsU(ix,ix +1) =   -imTime(1)*mu*16/(12*hy^2);
        lhsU(ix,ix -1) =   -imTime(1)*mu*16/(12*hy^2);
        lhsU(ix,ix +2) =   -imTime(1)*mu*(-1/(12*hy^2));
        lhsU(ix,ix -2) =   -imTime(1)*mu*(-1/(12*hy^2));
        
        lhsU(ix,ix + Nyg)   =   -imTime(1)*mu*16/(12*hx^2);
        lhsU(ix,ix - Nyg)   =   -imTime(1)*mu*16/(12*hx^2);
        lhsU(ix,ix + 2*Nyg) =   -imTime(1)*mu*(-1/(12*hx^2));
        lhsU(ix,ix - 2*Nyg) =   -imTime(1)*mu*(-1/(12*hx^2));
        
    end
    
end

for ix = M-3*Nyg + 1 : M-2*Nyg
    
    if BC == 2 || BC == 6
    
        lhsU(ix,ix) = 1;
        lhsU(ix,ix - M + 5*Nyg) = -1;
        
    else
        
    lhsU(ix,ix) = 1;
    
    end
end

for ix = M-2*Nyg + 1 : M-2*Nyg + 3

    if BC == 2 || BC == 6
    
        lhsU(ix,ix) = 1;
        lhsU(ix,ix - M + 5*Nyg) = -1;
        
    else
        
    lhsU(ix,ix) = 1;
    
    end
    
end

for ix = M-2*Nyg + 4 : M-Nyg-3
    
    if BC == 1
        
        lhsU(ix,ix        ) =  ( 8)/(12*hx);
        lhsU(ix,ix +   Nyg) =  (-1)/(12*hx);
        lhsU(ix,ix - 2*Nyg) =  (-8)/(12*hx);
        lhsU(ix,ix - 3*Nyg) =  ( 1)/(12*hx);
    
    elseif BC == 2 || BC == 6
        
        lhsU(ix,ix) = 1;
        lhsU(ix,ix - M + 5*Nyg) = -1;
        
    elseif BC == 4
        
        lhsU(ix,ix) = 1;
    
    end
    
end

for ix = M-Nyg-2 : M - Nyg + 3
    
    if BC == 2 || BC == 6
    
        lhsU(ix,ix) = 1;
        lhsU(ix,ix - M + 5*Nyg) = -1;
        
    else
        
    lhsU(ix,ix) = 1;
    
    end
end


for ix = M - Nyg + 4: M-3
    
    if BC == 1
        
        lhsU(ix,ix        ) =  -1/(12*hx^2);
        lhsU(ix,ix -   Nyg) =  16/(12*hx^2);
        lhsU(ix,ix - 2*Nyg) = -30/(12*hx^2);
        lhsU(ix,ix - 3*Nyg) =  16/(12*hx^2);
        lhsU(ix,ix - 4*Nyg) =  -1/(12*hx^2);
                
    elseif BC == 2 || BC == 6
        
        lhsU(ix,ix) = 1;
        lhsU(ix,ix - M + 5*Nyg) = -1;
        
    elseif BC == 4
        lhsU(ix,ix) = 1;
    end
end

for ix = M-2:M
    
    if BC == 2 || BC == 6
    
        lhsU(ix,ix) = 1;
        lhsU(ix,ix - M + 5*Nyg) = -1;
        
    else
        
    lhsU(ix,ix) = 1;
    
    end
    
end

lhsU = sparse(lhsU);
end


function L = setLaplacian(L,hx,hy,mu,imTime,M,pts,lPts)
  
% get the coeffcients for the laplacian dxx + dyy in fourth order
  coeffC  =  (1-imTime(1)*mu*(-(30/(12*hy^2)+30/(12*hx^2)))) * ones(1,lPts);
  coeffY1 = -imTime(1)*mu*16/(12*hy^2)   * ones(1,lPts);
  coeffY2 = -imTime(1)*mu*(-1/(12*hy^2)) * ones(1,lPts);
  coeffX1 = -imTime(1)*mu*16/(12*hx^2)   * ones(1,lPts);
  coeffX2 = -imTime(1)*mu*(-1/(12*hx^2)) * ones(1,lPts);
% put them o the corresponding locations in the sparse matrx
  L = L + sparse(pts,pts,coeffC,M,M);
% dyy
  L = L + sparse(pts,pts+1,coeffY1,M,M);
  L = L + sparse(pts,pts-1,coeffY1,M,M);
  L = L + sparse(pts,pts+2,coeffY2,M,M);
  L = L + sparse(pts,pts-2,coeffY2,M,M);
% dxx
  L = L + sparse(pts,pts+Nyg,coeffX1,M,M);
  L = L + sparse(pts,pts-Nyg,coeffX1,M,M);
  L = L + sparse(pts,pts+2*Nyg,coeffX2,M,M);
  L = L + sparse(pts,pts-2*Nyg,coeffX2,M,M);

end

function L = setDx(L,h,side,Stride,M,pts,lPts)

  coeffX1 = ones(1,length(lPts)) * ( 8)/(12*h) * (-1)^(side+1);
  coeffX2 = ones(1,length(lPts)) * ( 1)/(12*h) * (-1)^(side+1);
  
  shift1  =  ((-1)^side)*2*Stride; 
  shift2  =  ((-1)^side)*3*Stride;
  shift3  = -((-1)^side)*  Stride;
  
  L   = L + sparse(pts,pts         , coeffX1,M,M);
  L   = L + sparse(pts,pts + shift1,-coeffX1,M,M);
  L   = L + sparse(pts,pts + shift2, coeffX2,M,M);
  L   = L + sparse(pts,pts + shift3,-coeffX2,M,M);
  
end

function setDxx
  
end


function setDxxx(L,h,side,Stride,M,pts,lPts)
  
  coeffX1 = ones(1,length(lPts)) * 1/(2*h^3) * (-1)^(side+1);
  coeffX2 = ones(1,length(lPts)) * 1/(  h^3) * (-1)^(side+1);
  
  shift1  =  ((-1)^side)*  Stride;
  shift2  =  ((-1)^side)*3*Stride;
  shift3  =  ((-1)^side)*4*Stride;
  
  L   = L + sparse(pts,pts         , coeffX1,M,M);
  L   = L + sparse(pts,pts + shift1,-coeffX2,M,M);
  L   = L + sparse(pts,pts + shift2, coeffX2,M,M);
  L   = L + sparse(pts,pts + shift3,-coeffX1,M,M);

end


function L = setOne(L,M,pm,pts,lPts)

  coeff = pm*ones(1,length(lPts));
  L     = L + sparse(pts,pts,coeff,M,M);

end


