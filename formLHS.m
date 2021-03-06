function [lhsU,lhsV] = formLHS(fS)
%--------------------------------------------------------------------------
% sub functions contained in this code:
% 1. function L = setLaplacian(L,hx,hy,Nyg,mu,imTime,M,pts,lPts)
% 2. function L = setDx(L,h,side,Stride,M,pts,lPts)
% 3. function L = setDxx(L,h,side,Stride,M,pts,lPts)
% 4. function L = setDxxx(L,h,side,Stride,M,pts,lPts)
% 5. function L = setOne(L,M,pm,pts,matchPts,lPts)
% 6. function L = setZero(L,M,pts,lPts)
% 7. function L = setExt(L,side,Stride,M,pts)
%--------------------------------------------------------------------------

% get info from INS class

% if 1, use 4/6th order extrapolation on the first ghost line
% for tangetial velocity ( tangentExtOrder determines the order)
% 042317 fm
tangentExt      = fS.tangentExt;
tangentExtOrder = fS.tangentExtOrder;

Nxg = fS.Nxg;
Nyg = fS.Nyg;

Nx = fS.Nx;
Ny = fS.Ny;

ia = fS.ia;
ib = fS.ib;
ja = fS.ja;
jb = fS.jb;

hx  = fS.hx;
hy  = fS.hy;

BC  = fS.BC;

%extOrder   = fS.extOrder;
imTime     = fS.cBI;
mu         = fS.mu;

% total number of grid points(including ghosts)
M    = Nyg*Nxg;
% allocate LHS matrix for U
lhsU = spalloc(M,M,7*M);
lhsV = spalloc(M,M,7*M);

% get the indecies for all the interior points for U
[intPts,lIntPts] = getIndex(fS,ia+1,ib-1,ja+1,jb-1);
% put the laplacian operator onto these indecies
lhsU = setLaplacian(lhsU,hx,hy,Nyg,mu,imTime,M,intPts,lIntPts);
lhsV = setLaplacian(lhsV,hx,hy,Nyg,mu,imTime,M,intPts,lIntPts);

% assembling boundary conditions
% x-direction
for side = 0:1
    for axis = 0:1
        
        localBC = BC(axis+1,side+1);
        % on the boundary
        pos = 0;
        [bcPts,lBcPts] = getBCGLIndex(fS,axis,side,pos);
        
        bcxStart = (side==0)*ia + (side==1)*ib;
        bcyStart = (side==0)*ja + (side==1)*jb;
        
        bcxEnd = (axis==1)*((side==1)*ia + (side==0)*ib) +  (axis==0)*(bcxStart);
        bcyEnd = (axis==0)*((side==1)*ja + (side==0)*jb) +  (axis==1)*(bcyStart);
        
        if bcxStart>bcxEnd
            temp = bcxEnd;
            bcxEnd = bcxStart;
            bcxStart = temp;
        end
        
        if bcyStart>bcyEnd
            temp = bcyEnd;
            bcyEnd = bcyStart;
            bcyStart = temp;
        end
        
        [allBcPts,lAllBcPts] = getIndex(fS,bcxStart,bcxEnd,bcyStart,bcyEnd);
        lhsU = setZero(lhsU,M,allBcPts,lAllBcPts);
        lhsV = setZero(lhsV,M,allBcPts,lAllBcPts);
        
        switch localBC
            
            case 1
                % u = u(t)
                pm    = 1;
                lhsU  = setOne(lhsU,M,pm,bcPts,bcPts,lBcPts);
                lhsV  = setOne(lhsV,M,pm,bcPts,bcPts,lBcPts);
                
            case 2
                % Laplace(U) = ...
                if side == 0
                    lhsU  = setLaplacian(lhsU,hx,hy,Nyg,mu,imTime,M,bcPts,lBcPts);
                    lhsV  = setLaplacian(lhsV,hx,hy,Nyg,mu,imTime,M,bcPts,lBcPts);
                    
                elseif side == 1
                    
                    %[bcPts,lBcPts] = getBCGLIndex(fS,axis,side,pos);
                    
                    lhsU  = setOne(lhsU,M,+1,bcPts,bcPts,lBcPts);
                    lhsV  = setOne(lhsV,M,+1,bcPts,bcPts,lBcPts);
                    
                    matchSide = 0;
                    [matchPts,lMatchPts] =  getBCGLIndex(fS,axis,matchSide,-pos);
                    
                    lhsU  = setOne(lhsU,M,-1,bcPts,matchPts,lMatchPts);
                    lhsV  = setOne(lhsV,M,-1,bcPts,matchPts,lMatchPts);
                    
                end
                
            case 3
                
                % u = u(t)
                pm    = 1;
                lhsU  = setOne(lhsU,M,pm,bcPts,bcPts,lBcPts);
                lhsV  = setOne(lhsV,M,pm,bcPts,bcPts,lBcPts);
                
            case 4
                
            case 5
                
        end
        %----------------------------------------------------------------
        % on the ghost lines
        %----------------------------------------------------------------
        
        switch localBC
            
            case 1
                
                % u_x = -v_y
                
                pos = 1;
                
                [bcPts,lBcPts] = getBCGLIndex(fS,axis,side,pos);
                
                if axis == 0
                    lhsU  = setDx(lhsU,hx,side,Nyg,M,bcPts,lBcPts);
                    
                    if tangentExt == 1
                        
                        lhsV  = setExt(lhsV,side,Nyg,M,bcPts);
                        
                    elseif tangentExt == 0
                        
                        lhsV  = setDxx(lhsV,hx,side,Nyg,M,bcPts,lBcPts,pos);
                        
                    end
                    
                elseif axis == 1
                    
                    if tangentExt == 1
                        lhsU  = setExt(lhsU,side,1,M,bcPts);
                        
                    elseif tangentExt == 0
                        
                        lhsU  = setDxx(lhsU,hy,side,1,M,bcPts,lBcPts,pos);
                        
                    end
                    
                    lhsV  = setDx(lhsV,hy,side,1,M,bcPts,lBcPts);
                    
                end
                
                pos = 2;
                
                [bcPts,lBcPts] = getBCGLIndex(fS,axis,side,pos);
                
                if axis == 0
                    lhsU  = setDxxx(lhsU,hx,side,Nyg,M,bcPts,lBcPts);
                    lhsV  = setExt(lhsV,side,Nyg,M,bcPts);
                    
                elseif axis == 1
                    lhsU  = setExt(lhsU,side,1,M,bcPts);
                    lhsV  = setDxx(lhsV,hy,side,1,M,bcPts,lBcPts,pos);
                    
                end
                
            case 2
                % u(ia-i,j)  = u(Nx - i,j)
                
                for pos = 1:2;
                    
                    [bcPts,lBcPts] = getBCGLIndex(fS,axis,side,pos);
                    
                    lhsU  = setOne(lhsU,M,+1,bcPts,bcPts,lBcPts);
                    lhsV  = setOne(lhsV,M,+1,bcPts,bcPts,lBcPts);
                    
                    matchSide = abs(side-1);
                    [matchPts,lMatchPts] =  getBCGLIndex(fS,axis,matchSide,-pos);
                    
                    lhsU  = setOne(lhsU,M,-1,bcPts,matchPts,lMatchPts);
                    lhsV  = setOne(lhsV,M,-1,bcPts,matchPts,lMatchPts);
                    
                end
                
            case 3
                % u = u(t)
                for pos = 1:2;
                    
                    [bcPts,lBcPts] = getBCGLIndex(fS,axis,side,pos);
                    
                    pm    = 1;
                    lhsU  = setOne(lhsU,M,pm,bcPts,bcPts,lBcPts);
                    lhsV  = setOne(lhsV,M,pm,bcPts,bcPts,lBcPts);
                    
                end
                
            case 4
                
            case 5
                
        end
        
    end
end

% Note the righ-bottom corner point (A)
% _____A
%|     |
%|     |
%|_____|
%
% since the boundary condition is given in order of:
%left - bottom - right - top
% point A is contaminated by the top boundary condition

%----------------------------------------------------------------
% on the physical corner and its extension (ghost points)
%----------------------------------------------------------------

for sideX = 0:1
    for sideY = 0:1
        
        pxC = (sideX==0)*ia + (sideX==1)*ib;
        pyC = (sideY==0)*ja + (sideY==1)*jb;
        
        chooseBC(1) = BC(1,1+sideX);
        chooseBC(2) = BC(2,1+sideY);
        
        if max( chooseBC ) < 4
            
            [localBC,whichBC] =  max( chooseBC );
            axisBC = (whichBC -1);
            
        else
            fprintf('these could be wrong, since the developer havenot investigated this yet.\n');
            fprintf('proceed with extreme caution \n');
            pause;
            
        end
        
        switch localBC
            
            case 1
                %do nothing
            case 2
                
                side = (axisBC==0)*sideX + (axisBC==1)*sideY;
                sideMatch = abs(side-1);
                
                if side == 0
                    %do nothing
                    
                elseif side == 1
                    
                    % * = (Nx,Ny)
                    %
                    %       #
                    %       #
                    % xxxxxx*oo
                    %       x
                    %       x
                    %       x
                    
                    for pos = 0:2 % fixing *oo
                        
                        px = pxC + (axisBC == 0) * pos;
                        py = pyC + (axisBC == 1) * pos;
                        
                        pxIndex = getIndex(fS,px,px,py,py);
                        lPxIndex = 1;
                        lhsU = setZero(lhsU,M,pxIndex,lPxIndex);
                        lhsV = setZero(lhsV,M,pxIndex,lPxIndex);
                        
                        lhsU = setOne(lhsU,M,+1,pxIndex,pxIndex,lPxIndex);
                        lhsV = setOne(lhsV,M,+1,pxIndex,pxIndex,lPxIndex);
                        
                        matchPx = px - (axisBC==0)*(Nx-1);
                        matchPy = py - (axisBC==1)*(Ny-1);
                        matchPxIndex = getIndex(fS,matchPx,matchPx,matchPy,matchPy);
                        
                        lhsU = setOne(lhsU,M,-1,pxIndex,matchPxIndex,lPxIndex);
                        lhsV = setOne(lhsV,M,-1,pxIndex,matchPxIndex,lPxIndex);
                        
                    end
                    
                    
                    for pos = 1:2 % fixing ##
                        
                        px = pxC + (axisBC == 1) * (-1)^(sideX+1)* pos;
                        py = pyC + (axisBC == 0) * (-1)^(sideY+1)* pos;
                        
                        pxIndex = getIndex(fS,px,px,py,py);
                        lPxIndex = 1;
                        
                        lhsU = setZero(lhsU,M,pxIndex,lPxIndex);
                        lhsV = setZero(lhsV,M,pxIndex,lPxIndex);
                        
                        lhsU = setOne(lhsU,M,+1,pxIndex,pxIndex,lPxIndex);
                        lhsV = setOne(lhsV,M,+1,pxIndex,pxIndex,lPxIndex);
                        
                        matchPx = px - (axisBC==0)*(Nx-1);
                        matchPy = py - (axisBC==1)*(Ny-1);
                        matchPxIndex = getIndex(fS,matchPx,matchPx,matchPy,matchPy);
                        
                        lhsU = setOne(lhsU,M,-1,pxIndex,matchPxIndex,lPxIndex);
                        lhsV = setOne(lhsV,M,-1,pxIndex,matchPxIndex,lPxIndex);
                        
                    end
                    
                    
                end
                
                
            case 3
                %do nothing
                
        end
        
        
    end
end


%----------------------------------------------------------------
% now fix em ghost corner points
%----------------------------------------------------------------
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
        
        px1Index = getIndex(fS,px1,px1,py1,py1);
        px2Index = getIndex(fS,px2,px2,py2,py2);
        px3Index = getIndex(fS,px3,px3,py3,py3);
        px4Index = getIndex(fS,px4,px4,py4,py4);
        
        chooseBC(1) = BC(1,1+sideX);
        chooseBC(2) = BC(2,1+sideY);
        
        if max( chooseBC ) < 4
            
            [localBC,whichBC] =  max( chooseBC );
            axisBC = (whichBC -1);
            
        else
            fprintf('these could be wrong, since the developer havenot investigated this yet.\n');
            fprintf('proceed with extreme caution \n');
            pause;
            
        end
        
        switch localBC
            
            case 1
                
                for i = 1:3
                    px1E(i) = px1 + (-1)^(sideX)*i;
                    py1E(i) = py1 + (-1)^(sideY)*i;
                    
                    px2E(i) = px2 + (-1)^(sideX)*i*(2);
                    py2E(i) = py2 + (-1)^(sideY)*i;
                    
                    px3E(i) = px3 + (-1)^(sideX)*i;
                    py3E(i) = py3 + (-1)^(sideY)*i*(2);
                    
                    px4E(i) =  px1E(i);
                    py4E(i) =  py1E(i);
                end
                
                coeff  =[1,-15/4 ,3 ,-1/4];
                coeff4 =[1,-30   ,32,-3  ];
                
                %                 px1E=[];
                %                 px2E=[];
                %                 px3E=[];
                %                 px4E=[];
                %
                %                 py1E=[];
                %                 py2E=[];
                %                 py3E=[];
                %                 py4E=[];
                %
                %                 coeff=1;
                %                coeff4=1;
                
                
            case 2
                
                stripeMatch = Nx*(axisBC==0) + Ny*(axisBC==1);
                
                side = sideX*(axisBC==0) + sideY*(axisBC==1);
                
                px1E = px1 + (-1)^side*(stripeMatch - 1);
                py1E = py1 ;
                
                px2E = px2 + (-1)^side*(stripeMatch - 1);
                py2E = py2 ;
                
                px3E = px3 + (-1)^side*(stripeMatch - 1);
                py3E = py3 ;
                
                px4E = px4 + (-1)^side*(stripeMatch - 1);
                py4E = py4 ;
                
                coeff =[1,-1];
                coeff4=[1,-1];
                
                
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
                coeff4=1;
                
                
                
        end
        
        px1 = [px1,px1E];
        px2 = [px2,px2E];
        px3 = [px3,px3E];
        px4 = [px4,px4E];
        
        py1 = [py1,py1E];
        py2 = [py2,py2E];
        py3 = [py3,py3E];
        py4 = [py4,py4E];
        
        for i = 1:length(px1)
            p1Index = getIndex(fS,px1(i),px1(i),py1(i),py1(i));
            p2Index = getIndex(fS,px2(i),px2(i),py2(i),py2(i));
            p3Index = getIndex(fS,px3(i),px3(i),py3(i),py3(i));
            p4Index = getIndex(fS,px4(i),px4(i),py4(i),py4(i));
            
            lhsU    = lhsU + sparse(px1Index,p1Index,coeff(i),M,M);
            lhsU    = lhsU + sparse(px2Index,p2Index,coeff(i),M,M);
            lhsU    = lhsU + sparse(px3Index,p3Index,coeff(i),M,M);
            lhsU    = lhsU + sparse(px4Index,p4Index,coeff4(i),M,M);
            
            lhsV    = lhsV + sparse(px1Index,p1Index,coeff(i),M,M);
            lhsV    = lhsV + sparse(px2Index,p2Index,coeff(i),M,M);
            lhsV    = lhsV + sparse(px3Index,p3Index,coeff(i),M,M);
            lhsV    = lhsV + sparse(px4Index,p4Index,coeff4(i),M,M);
        end
        
    end
end


end


function L = setLaplacian(L,hx,hy,Nyg,mu,imTime,M,pts,lPts)

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

coeffX1 = ones(1,lPts) * ( 8)/(12*h) * (-1)^(side+1);
coeffX2 = ones(1,lPts) * ( 1)/(12*h) * (-1)^(side+1);

shift1  =  ((-1)^side)*2*Stride;
shift2  =  ((-1)^side)*3*Stride;
shift3  = -((-1)^side)*  Stride;

L   = L + sparse(pts,pts         , coeffX1,M,M);
L   = L + sparse(pts,pts + shift1,-coeffX1,M,M);
L   = L + sparse(pts,pts + shift2, coeffX2,M,M);
L   = L + sparse(pts,pts + shift3,-coeffX2,M,M);

end

function L = setDxx(L,h,side,Stride,M,pts,lPts,pos)

coeffX1 = ones(1,lPts) * (  16)/(12*h^2) ;
coeffX2 = ones(1,lPts) * ( -30)/(12*h^2) ;
coeffX3 = ones(1,lPts) * ( - 1)/(12*h^2) ;

awayFromBC = pos-2;

shift1  =  ((-1)^side)*(awayFromBC  )*Stride;
shift2  =  ((-1)^side)*(awayFromBC+1)*Stride;
shift3  =  ((-1)^side)*(awayFromBC+2)*Stride;
shift4  =  ((-1)^side)*(awayFromBC+3)*Stride;
shift5  =  ((-1)^side)*(awayFromBC+4)*Stride;

%L   = L + sparse(pts,pts         , coeffX1,M,M);
L   = L + sparse(pts,pts + shift1, coeffX3,M,M);
L   = L + sparse(pts,pts + shift2, coeffX1,M,M);
L   = L + sparse(pts,pts + shift3, coeffX2,M,M);
L   = L + sparse(pts,pts + shift4, coeffX1,M,M);
L   = L + sparse(pts,pts + shift5, coeffX3,M,M);


end


function L = setDxxx(L,h,side,Stride,M,pts,lPts)

coeffX1 = ones(1,lPts) * 1/(2*h^3) * (-1)^(side+1);
coeffX2 = ones(1,lPts) * 1/(  h^3) * (-1)^(side+1);

shift1  =  ((-1)^side)*  Stride;
shift2  =  ((-1)^side)*3*Stride;
shift3  =  ((-1)^side)*4*Stride;

L   = L + sparse(pts,pts         , coeffX1,M,M);
L   = L + sparse(pts,pts + shift1,-coeffX2,M,M);
L   = L + sparse(pts,pts + shift2, coeffX2,M,M);
L   = L + sparse(pts,pts + shift3,-coeffX1,M,M);

end


function L = setOne(L,M,pm,pts,matchPts,lPts)

coeff = pm*ones(1,lPts);
L     = L + sparse(pts,matchPts,coeff,M,M);

end

function L = setZero(L,M,pts,lPts)

coeffZero = zeros(lPts,M);

L(pts,:)     = coeffZero;

end

function L = setExt(L,side,Stride,M,pts)

extOrder = 4;
%extOrder = fS.tangentOrder;
if extOrder == 6
    coeff = [1, -6, 15, -20, 15, -6, 1];
elseif extOrder == 4
    coeff = [1, -4, 6, -4, 1];
end
    


for i = 1:extOrder+1
    
    shift  =  ((-1)^side)*(i-1)*Stride;
    
    ptsShift = pts + shift;
    
    L   = L + sparse(pts,ptsShift, coeff(i),M,M);
end

end
