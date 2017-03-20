function miscellaneous

% calculate the coefficient for 2nd order extrapolation in time 
%--------------------------------------------------------------
syms t2 t1 tp1 tp2

A = [1 1 1;...
    t1 tp1 tp2;...
    t1^2 tp1^2 tp2^2];

b = [1;t2;t2^2];

x  = A\b;
%--------------------------------------------------------------

% check if that returns 3 -3 1 given equidistant dt
%--------------------------------------------------------------
tp2 = 1;
tp1 = 2;
t1  = 3;
t2 = 4;

x(1)=         -(t2*tp1 + t2*tp2 - tp1*tp2 - t2^2)/((t1 - tp1)*(t1 - tp2));
x(2)=           (t1*t2 - t1*tp2 + t2*tp2 - t2^2)/((t1 - tp1)*(tp1 - tp2));
x(3)=-(t1*t2 - t1*tp1 + t2*tp1 - t2^2)/(t1*tp1 - t1*tp2 - tp1*tp2 + tp2^2);
%--------------------------------------------------------------




end