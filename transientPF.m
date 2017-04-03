mu = 0.1;
p0 = 1;
y = linspace(-1,1,100);
h = 1;
k = 1;

p1 = 1;
t = 0;

U0 = p0/(2*mu)*(h^2-y.^2);

U1 = 1i*p1/(k)*(cosh(y*sqrt(1i*k/mu))/cosh(h*sqrt(1i*k/mu)) - 1)*exp(1i*k*t);


plot(real(U1),y)
