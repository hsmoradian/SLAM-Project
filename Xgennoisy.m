function xhist = Xgennoisy(F,G,Gamma,r,Q,xhat0,kmax,N)
x(:,1) = xhat0;
for k = 1:kmax
    vrnd = randn(N,1);
    v(:,k) = chol(Q)*vrnd;
    x(:,k+1) = F*x(:,k) + G*r(:,k) + Gamma*v(:,k);
end
xhist = x(:,2:end);