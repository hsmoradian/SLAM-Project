function e = egen(Xkftilde,Pkf,kmax)
for k = 1:kmax
    N = 5;
    [S,P] = qr(Pkf(:,:,k));
    e(1,k) = (Xkftilde(:,k)'*S)*(P'*Xkftilde(:,k));
end