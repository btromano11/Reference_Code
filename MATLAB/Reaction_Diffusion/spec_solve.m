function solfvecsol = spec_solve()

% parameters
tspan = 0:0.5:4;
beta = 1;
D1 = .1;
D2 = .1;
m = 1;
L=20;
n=64;

%domain
x2 = linspace(-L/2,L/2,n+1);x=x2(1:n);
y=x;
[X,Y]=meshgrid(x,y);

%spectral variables
kx=(2*pi/L)*[0:(n/2-1) (-n/2):-1];kx(1)=10^-6;
ky=kx;
[KX,KY]=meshgrid(kx,ky);
K = KX.^2+KY.^2;

%initial conditions
uinit = tanh(sqrt(X.^2+Y.^2)).*cos(m*angle(X+1i*Y)-(sqrt(X.^2+Y.^2)));
vinit = tanh(sqrt(X.^2+Y.^2)).*sin(m*angle(X+1i*Y)-(sqrt(X.^2+Y.^2)));

ufinit = fft2(uinit);
ufvecinit = reshape(ufinit,n^2,1);

vfinit = fft2(vinit);
vfvecinit = reshape(vfinit,n^2,1);

solfvecinit = [ufvecinit; vfvecinit];

[t,solfvecsol]=ode45(@(t,solfvec) spec_rhs(t,solfvec,beta,D1,D2,K,n),tspan,solfvecinit);


% for j=1:length(t)
%     
%     wfvecsol = solfvecsol(j,1:n^2);
%     curw=real(ifft2(reshape(wfvecsol(1,:),n,n)));
%     pcolor(X,Y,curw);shading interp;
%     drawnow;
%     
%     pause(0.2);
% end

