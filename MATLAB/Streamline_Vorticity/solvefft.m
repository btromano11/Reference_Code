function w = solvefft(L,n)

nsq = n^2;

%create x/y grid, reshape to vector
x2=linspace(-L,L,n+1);
x=x2(1:n);
y2=linspace(-L,L,n+1);
y=y2(1:n);
[X, Y] = meshgrid(x,y);
winit = exp(-X.^2-(Y.^2)/20);
w = winit;

%allocate phi
phi = 1:nsq;
phi = phi';

%fft parameters
kx=(pi/(L))*[0:(n/2-1) (-n/2):-1];kx(1)=10^-6;
ky=kx;
[X,Y]=meshgrid(x,y);
[KX,KY]=meshgrid(kx,ky);
Kvec = KX.^2+KY.^2;

%calculate constants
delta = x(2)-x(1);
tspan = 0:0.5:4;
v = 0.001;

%create matrix operators
A = dx2dy2(n, delta);
A(1,1) = 2;
B = dx(n, delta);                                                                                                                                         
C = dy(n, delta);

%integrate in time
[t,w] = ode45(@(t,w) rhsfft(t,w,phi,A,B,C,v,Kvec,nsq),tspan,winit);


