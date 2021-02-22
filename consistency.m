function e = consistency(alpha,N,Xkftilde,Pkf,kmax)
for k = 1:kmax
    e(1,k) = Xkftilde(:,k+1)'*Pkf(:,:,k)*Xkftilde(:,k+1);
end