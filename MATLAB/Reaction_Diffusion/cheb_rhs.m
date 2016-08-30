function rhs = cheb_rhs(t,solvec,beta,D1,D2,Lap,n,L)

u = solvec(1:(n+1)^2);
v = solvec((n+1)^2+1:2*(n+1)^2);

A2 = u.^2+v.^2;

unew = (1-A2).*u - (-beta*A2.*v) + D1*Lap*u;
vnew =  (-beta*A2.*u) + (1-A2).*v + D2*Lap*v;

rhs = [unew;vnew];