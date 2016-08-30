function A = dx2dy2(n, delta)

nsq = n^2;

%Fill diagonals
for i = 1:n:nsq
    B(1,i:i+(n-1)) = ones(n,1);
    B(2,i:i+(n-1)) = ones(n,1);
    B(3, i) = 1;
    B(3,i+1:i+n-1) = zeros(n-1,1);
    B(4, i:i+n-2) = (ones(n-1,1));
    B(4,i+n-1) = 0;
    B(5,i:i+(n-1)) = -4*(ones(n,1));
    B(6, i) = 0;
    B(6,i+1:i+n-1) = ones(n-1,1);
    B(7,i:i+n-2) = zeros(n-1,1);
    B(7, i+n-1) = 1;
    B(8,i:i+(n-1)) = ones(n,1);
    B(9,i:i+(n-1)) = ones(n,1);
end

%multiply by constant, make diagonal column vectors
B = B';
B = B/(delta^2);

%create sparse matrix
A = spdiags(B,[-(nsq-n) -n -(n-1) -1 0 1 (n-1) n (nsq-n)],nsq,nsq);



