function zhist = Zgen(xnoisy,R,kmax,N)
for k = 1:kmax
    wrnd = randn(N,1);
    w(:,k) = chol(R)*wrnd;
    z(:,k) = xnoisy(N+1:2*N,k) + w(:,k);
end
zhist = z;