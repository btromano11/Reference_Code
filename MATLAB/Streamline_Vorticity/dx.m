function A = dx(n, delta)

nsq = n^2;

%create diagonals, multiply by constant, make diagonal column vectors
B = ones(nsq,1);

B = B/(2*delta);

%create sparse matrix
A = spdiags([B -B B -B],[-(nsq-n) -n n (nsq-n)],nsq,nsq);
