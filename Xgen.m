function xhist = Xgen(F,G,r,xhat0,kmax,N)
x(:,1) = xhat0;
for k = 1:kmax
    x(:,k+1) = F*x(:,k) + G*r(:,k);
end
xhist = x;