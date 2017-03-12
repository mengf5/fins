function [lhsU] = formLHS(fS)

% get info from INS class 
Nxg = fS.Nxg;
Nyg = fS.Nyg;

Nx = fS.Nx;
Ny = fS.Ny;


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
    
    
    switch localBC
	       
	  case 1

	    for i = 1:3
	      px1E(i) = px1 + (-1)^(sideX)*i
	      py1E(i) = py1 + (-1)^(sideY)*i

	      px2E(i) = px1 + (-1)^(sideX)*i*(2)
	      py2E(i) = py1 + (-1)^(sideY)*i

	      px3E(i) = px1 + (-1)^(sideX)*i
	      py3E(i) = py1 + (-1)^(sideY)*i*(2)

	      px4E(i) = px1 + (-1)^(sideX)*i*(2)
	      py4E(i) = py1 + (-1)^(sideY)*i*(2)
	    end
	    
	    coeff=[1,-15/4,3,-1/4];


	  case 2

	    px1E = px1 + (-1)^sideX*(Nx - 1);
	    py1E = py1 + (-1)^sideX*(Ny - 1);

	    px2E = px2 + (-1)^sideX*(Nx - 1);
	    py2E = py2 + (-1)^sideX*(Ny - 1);

	    px3E = px3 + (-1)^sideX*(Nx - 1);
	    py3E = py3 + (-1)^sideX*(Ny - 1);

	    px4E = px4 + (-1)^sideX*(Nx - 1);
	    py4E = py4 + (-1)^sideX*(Ny - 1);

	    coeff=[1,-1];


	  case 3

	    px1E=[];
	    px2E=[];
	    px3E=[];
	    px4E=[];

	    py1E=[];
	    py2E=[];
	    py3E=[];
	    py4E=[];

	    coeff=1;



	end

	px1 = [px1;px1E];
	px2 = [px2;px2E];
	px3 = [px3;px3E];
	px4 = [px4;px4E];

	py1 = [py1;py1E];
	py2 = [py2;py2E];
	py3 = [py3;py3E];
	py4 = [py4;py4E];
	
	p1Index = getIndex(fS,px1,px1,py1,py1);
	p2Index = getIndex(fS,px1,px1,py1,py1);
	p3Index = getIndex(fS,px1,px1,py1,py1);
	p4Index = getIndex(fS,px1,px1,py1,py1);
	
	L    = L + sparse(p1Index,p1Index,coeff,M,M);
	L    = L + sparse(p2Index,p2Index,coeff,M,M);
	L    = L + sparse(p3Index,p3Index,coeff,M,M);
	L    = L + sparse(p4Index,p4Index,coeff,M,M);
	
      end
    end    
    
    
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

function setDxx(L,h,side,Stride,M,pts,lPts)

  
  lhsU(ix,ix - 1) =   -1/(12*hy^2);
  lhsU(ix,ix    ) =   16/(12*hy^2);
  lhsU(ix,ix + 1) =  -30/(12*hy^2);
  lhsU(ix,ix + 2) =   16/(12*hy^2);
  lhsU(ix,ix + 3) =   -1/(12*hy^2);
	    
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

function L  = setExt(L,h,side,Stride,M,pts,lPts)

  extOrder = 6;
  coeff = [1, -6, 15, -20, 15, -6, 1];
  
  for i = 0:extOrder
    shift(i)  =  ((-1)^side)*i*Stride;
    pts       = pts
  end
  
  pts      = ones(1,1+extOrder) * pts;
  ptsShift = pts + shift;
  
  L   = L + sparse(pts,ptsShift, ptsShift,M,M);
  
end
