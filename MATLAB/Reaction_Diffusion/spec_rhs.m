function rhs = spec_rhs(t,solfvec,beta,D1,D2,K,n)

%reshape u,v have fourier and real space versions
fu = reshape(solfvec(1:n^2),[n,n]);
fv = reshape(solfvec(n^2+1:2*n^2),[n,n]);
u = real(ifft2(fu));
v = real(ifft2(fv));

lapu = -K.*fu;
lapv = -K.*fv;

A2 = u.^2+v.^2;

%calculate new u,v in real space
unew =  (1-A2).*u - (-beta*A2.*v) + D1*real(ifft2(lapu));
vnew =  (-beta*A2.*u) + (1-A2).*v + D2*real(ifft2(lapv));

rhs = [reshape(fft2(unew),n^2,1);reshape(fft2(vnew),n^2,1)];