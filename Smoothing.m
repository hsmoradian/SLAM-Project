function [xshat,Ps,xfhat,Pfbar] = Smoothing(F,Gamma,H,Q,R,xhat0,P0,thist,zhist,N)
xhat = xhat0;
P = P0;
for k = 1:1:N
    %Prediction
    xbar = F*xhat;
    Pbar = F*P*F'+Gamma*Q*Gamma';
    %Update
    v = zhist(k)-H*xbar;
    s = H*Pbar*H'+R;
    w = Pbar*H'/s;
    xhat = xbar+w*v;
    P = Pbar-w*s*w';
    
    xfhat(:,k) = xhat;
    xfbar(:,k) = xbar;
    Pfbar(:,:,k) = Pbar;
    Pf(:,:,k) = P;
end
xshat(:,N) = xfhat(:,N);
Ps(:,:,N) = Pf(:,:,N);
for k = N-1:-1:1
    C = Pf(:,:,k)*F'*inv(Pfbar(:,:,k));
    xshat(:,k) = xfhat(:,k) + C*(xshat(:,k+1)-xfbar(:,k+1));
    Ps(:,:,k) = Pf(:,:,k)+C*(Ps(:,:,k+1)-Pfbar(:,:,k))*C';
end