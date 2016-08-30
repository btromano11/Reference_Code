function A = dy(n, delta)

nsq = n^2;

%generate diagonals
for i = 1:n:nsq
    B(1, i) = 1;
    B(1,i+1:i+n-1) = zeros(n-1,1);
    B(2, i:i+n-2) = -1*(ones(n-1,1));
    B(2,i+n-1) = 0;
    B(3, i) = 0;
    B(3,i+1:i+n-1) = ones(n-1,1);
    B(4,i:i+n-2) = zeros(n-1,1);
    B(4, i+n-1) = -1;
end

%multiply by constant, make diagonal column vectors
B = B';
B = B/(2*delta);

%make sparse matrix
A = spdiags([B],[-(n-1) -1 1 (n-1)],nsq,nsq);




