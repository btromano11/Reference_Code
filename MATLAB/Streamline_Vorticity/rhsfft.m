function rhs=rhsfft(t,w,phi,A,B,C,v,Kvec,nsq,n)

%solve elliptic problem
%reshape omega
w = reshape(w,[n,n]);

%transform to fourier, solve for phi
wfft = fft2(w);
phifft = wfft./(-Kvec);
phi = real(ifft2(phifft));

%reshape phi and omega
phi = reshape(phi,nsq,1);
w = reshape(w,nsq,1);

%evaluate rhs
rhs= (-B*phi).*(C*w) + (C*phi).*(B*w) + v*(A*w);