function [Xkf,Pkf,e] = KalmanFilterloc(F,r,sigmav2,sigmaw2,xhat0,P0loc,zhist,kmax,N,d,Aj,a,b,dT,xnoisy)
xhat = xhat0;
for i = 1:N
    P(:,:,i) = P0loc;
end
for k = 1:kmax
    for i = 1:N
        Gamma = [-a*b*d(i);b*d(i)];
        %Prediction
        input = dT*[-a*b*Aj(i,:)*xnoisy(N+1:2*N,k);r(i,k)+b*Aj(i,:)*xnoisy(N+1:2*N,k)];
        xbari = F(:,:,i)*[xhat(i,1);xhat(i+N,1)]+input;
        xbar(i,1) = xbari(1,1);
        xbar(i+N,1) = xbari(2,1);
        Pbar(:,:,i) = F(:,:,i)*P(:,:,i)*F(:,:,i)'+Gamma*sigmav2*Gamma';
        %Update
        v = zhist(i,k)-xbar(i+N,1);
        s = [0 1]*Pbar(:,:,i)*[0;1]+sigmaw2;
        w = Pbar(:,:,i)*[0;1]/s;
        xhati = [xbar(i,1);xbar(i+N,1)]+w*v;
        xhat(i,1) = xhati(1,1);
        xhat(i+N,1) = xhati(2,1);
        P(:,:,i) = Pbar(:,:,i)-w*s*w';
        %Data
        Xkf(i,k) = xhat(i,1);
        Xkf(i+N,k) = xhat(i+N,1);
        Pkf(:,:,i,k) = P(:,:,i);
        Pbarkf(:,:,i,k) = Pbar(:,:,i);
        e(i,1,k) = v'*inv(s)*v;
    end
end