function fLCIMEXBDF(varargin)

% r_i = \f{t^{n+1-i}}{t^{n+1}}
% h = t^{n+1} - t^{n}
syms r1 r2 r3 r4 r5;
% undetermined coefficients

%
% \[ u^{n+1} + \sum_{i=1}^{i = s} a u^{n+1-i} = fLC(s)*h* \p_t u^{n+1} + b*h* \p_t u^{n}\]
%

order = 4;
localOrder = 5;
addTerm = 'L';
fLCBDF =0 ;

dosplayCoeff = 0;
if nargin == 1
   
    order = varargin{1};
   
end

if nargin == 2
   
    dosplayCoeff = varargin{2};
   
end


%BDF

if order == 2
    
    fLC = 2/3;
    
    if localOrder==3
        
        A = [1    1   0; ...
            r1   r2   -(1-r1); ...
            r1^2 r2^2 -2*r1*(1-r1)];
        
        B = [-1;
            fLC*(1-r1)-1;
            fLC*(1-r1)*2-1];
        
        x = linsolve(A, B);
    
    elseif localOrder == 2
        
        A = [1    1 ; ...
            r1   r2   ...
            ];
        
        B = [-1;
            fLC*(1-r1)-1
            ];
        
        x = linsolve(A, B);
        
    end
    
    
    
    
    
    %display(x);
    
elseif order == 4
    
    fLC = 12/25;
    if localOrder == 4
        
        A = [
            1    1    1     1     ;
            r1   r2    r3    r4   ;
            r1^2 r2^2  r3^2  r4^2 ;
            r1^3 r2^3  r3^3  r4^3
            ];
        
        B = [
            -1;
            fLC*(1-r1)-1;
            fLC*(1-r1)*2-1;
            fLC*(1-r1)*3-1
            ];
        
    elseif localOrder == 5
        
        if fLCBDF == 0
            
            A = [
                1    1    1     1             0 ;
                r1   r2    r3    r4     -(1-r1) ;
                r1^2 r2^2  r3^2  r4^2   -2*(1-r1) ;
                r1^3 r2^3  r3^3  r4^3   -3*(1-r1) ;
                r1^4 r2^4  r3^4  r4^4   -4*(1-r1)
                ];
            
            B = [
                -1;
                -1;
                -1;
                -1;
                -1
                ];
            
            
            
            
        elseif fLCBDF == 1
            
            if strcmp(addTerm,'R')
                A = [
                    1    1    1     1     0;
                    r1   r2    r3    r4   -(1-r1);
                    r1^2 r2^2  r3^2  r4^2 -2*r1*(1-r1);
                    r1^3 r2^3  r3^3  r4^3 -3*r1^2*(1-r1);
                    r1^4 r2^4  r3^4  r4^4 -4*r1^3*(1-r1);
                    ];
                
                B = [
                    -1;
                    fLC*(1-r1)-1;
                    fLC*(1-r1)*2-1;
                    fLC*(1-r1)*3-1;
                    fLC*(1-r1)*4-1
                    ];
                
            elseif strcmp(addTerm,'L')
                
                A = [
                    1    1    1     1     1;
                    r1   r2    r3    r4   r5;
                    r1^2 r2^2  r3^2  r4^2 r5^2;
                    r1^3 r2^3  r3^3  r4^3 r5^3;
                    r1^4 r2^4  r3^4  r4^4 r5^4
                    ];
                
                B = [
                    -1;
                    fLC*(1-r1)-1;
                    fLC*(1-r1)*2-1;
                    fLC*(1-r1)*3-1;
                    fLC*(1-r1)*4-1;
                    ];
                
                
            end
        end
    
    end
    
    x = linsolve(A, B);
    %display(x);
        
end

    checkConstDt(x);

%IMEXBDF

if order == 2
    
    A = [
        1-r1 1-r1; ...
        (1-r1)*2*r1   (1-r1)*2*r2; ...
        ];
    
    B = [
        1 + x(1)*r1 + x(2)*r2; 
        1 + x(1)*r1^2 + x(2)*r2^2;
        ];
    
    xE = linsolve(A, B);
    
    %display(xE);
    
elseif order == 4
    
    if fLCBDF == 0
        
        A = [
            1-r1 1-r1 1-r1 1-r1; ...
            (1-r1)*2*r1   (1-r1)*2*r2 (1-r1)*2*r3 (1-r1)*2*r4; ...
            (1-r1)*3*r1^2   (1-r1)*3*r2^2 (1-r1)*3*r3^2 (1-r1)*3*r4^2;
            (1-r1)*4*r1^3   (1-r1)*4*r2^3 (1-r1)*4*r3^3 (1-r1)*4*r4^3;
            ];
        
        B = [
            1 + x(1)*r1 + x(2)*r2 + x(3)*r3 + x(4)*r4;
            1 + x(1)*r1^2 + x(2)*r2^2 + x(3)*r3^2 + x(4)*r4^2;
            1 + x(1)*r1^3 + x(2)*r2^3 + x(3)*r3^3 + x(4)*r4^3;
            1 + x(1)*r1^4 + x(2)*r2^4 + x(3)*r3^4 + x(4)*r4^4;
            ];
        
        xE = linsolve(A, B);
        
        
        
    elseif fLCBDF == 1
        
        if strcmp(addTerm,'R')
            
            A = [
                1-r1 1-r1 1-r1 1-r1; ...
                (1-r1)*2*r1   (1-r1)*2*r2 (1-r1)*2*r3 (1-r1)*2*r4; ...
                (1-r1)*3*r1^2   (1-r1)*3*r2^2 (1-r1)*3*r3^2 (1-r1)*3*r4^2;
                (1-r1)*4*r1^3   (1-r1)*4*r2^3 (1-r1)*4*r3^3 (1-r1)*4*r4^3;
                ];
            
            B = [
                1 + x(1)*r1 + x(2)*r2 + x(3)*r3 + x(4)*r4;
                1 + x(1)*r1^2 + x(2)*r2^2 + x(3)*r3^2 + x(4)*r4^2;
                1 + x(1)*r1^3 + x(2)*r2^3 + x(3)*r3^3 + x(4)*r4^3;
                1 + x(1)*r1^4 + x(2)*r2^4 + x(3)*r3^4 + x(4)*r4^4;
                ];
            
            xE = linsolve(A, B);
            
        elseif strcmp(addTerm,'L')
            
            A = [
                1-r1 1-r1 1-r1 1-r1; ...
                (1-r1)*2*r1   (1-r1)*2*r2 (1-r1)*2*r3 (1-r1)*2*r4; ...
                (1-r1)*3*r1^2   (1-r1)*3*r2^2 (1-r1)*3*r3^2 (1-r1)*3*r4^2;
                (1-r1)*4*r1^3   (1-r1)*4*r2^3 (1-r1)*4*r3^3 (1-r1)*4*r4^3;
                ];
            
            B = [
                1 + x(1)*r1 + x(2)*r2 + x(3)*r3 + x(4)*r4 + x(5)*r5;
                1 + x(1)*r1^2 + x(2)*r2^2 + x(3)*r3^2 + x(4)*r4^2 + x(5)*r5^2;
                1 + x(1)*r1^3 + x(2)*r2^3 + x(3)*r3^3 + x(4)*r4^3 + x(5)*r5^3;
                1 + x(1)*r1^4 + x(2)*r2^4 + x(3)*r3^4 + x(4)*r4^4 + x(5)*r5^4;
                ];
            
            xE = linsolve(A, B);
            
        end
    end
    
end

checkConstDt(xE);


if dosplayCoeff == 1
    fprintf('BDF, alpha, beta are:: \n')
    display(x)
    
    fprintf('IMEX BDF, beta are:: \n')
    display(xE)
end
end


function checkConstDt(xI)

for i = 1:length(xI)
 
    f{i} = matlabFunction(xI(i));

end

if length(xI)<=3

    for i = 1:length(xI)
        f{i}(.9,.8)
    end
elseif length(xI)>3

    for i = 1:length(xI)
        f{i}(.9,.8,.7,.6)
    end
    
end    


end