function solvecsol = cheb_solve()

tspan = 0:0.5:4;
beta = 1;
D1 = .1;
D2 = .1;
L=20;
n=30;

%cheb matrix
[C,z] = cheb(n);
C2= C^2;
C2(1,:) = zeros(1,n+1);
C2(n+1,:) = zeros(1,n+1);

% create and rescale laplacian
I = eye(length(C2));
Lap = ((2/L)^2)*(kron(C2,I)+kron(I,C2));

%rescale and create domain
x=L*z/2;
y=x;
[X,Y]=meshgrid(x,y);

%initial conditions
m = 1;
uinit = tanh(sqrt(X.^2+Y.^2)).*cos(m*angle(X+1i*Y)-(sqrt(X.^2+Y.^2)));
uinitvec = reshape(uinit,(n+1)^2,1);
vinit = tanh(sqrt(X.^2+Y.^2)).*sin(m*angle(X+1i*Y)-(sqrt(X.^2+Y.^2)));
vinitvec = reshape(vinit,(n+1)^2,1);

%stack matrices
solvecinit = [uinitvec;vinitvec];

[t,solvecsol]=ode45(@(t,solvec) cheb_rhs(t,solvec,beta,D1,D2,Lap,n,L),tspan,solvecinit);

