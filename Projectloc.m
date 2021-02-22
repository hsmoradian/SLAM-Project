clc
clear all

N = 5;
a = 1;
b = 1;
I = eye(N);
PI = I - 1/N*ones(N,1)*ones(1,N);
dT = 0.01;
Aj = [0 1 0 0 1;
      1 0 1 0 0;
      0 1 0 1 0;
      0 0 1 0 1;
      1 0 0 1 0];
D = diag([2 2 2 2 2]);
L = D-Aj;
On = zeros(N,N);
for i = 1:N
    d(i) = sum(Aj(i,:));
    Aloc(:,:,i) = [0 a*b*d(i);-1 -a-b*d(i)];
    Floc(:,:,i) = eye(2)+dT*Aloc(:,:,i);
end
A = [On a*b*L;
     -I -a*I-b*L];
B = [On;a*I];
Sigmavm = [-a*b*Aj;b*Aj];
F = eye(2*N)+dT*A;
G = dT*B;
Gamma = dT*Sigmavm;
sigmav2 = .0001;
sigmaw2 = .01;
v = sigmav2*ones(1,N);
w = sigmaw2*ones(1,N);
Q = diag(v);
R = diag(w);
H = [On I];
intervalx0 = 3;
x0 = [zeros(N,1);intervalx0*randn(N,1)];
xhat0 = [zeros(N,1);intervalx0*randn(N,1)];
P0 = 10*eye(2*N);
P0loc = [0.00001 0;0 10];
tmax = 25;
kmax = tmax/dT;
for i = 1:kmax
    r(:,i) = rfun(i,dT);
    rbar(1,i) = mean(r(:,i));
end
Run = 5;
for run = 1:Run
    xnoisy = Xgennoisy(F,G,Gamma,r,Q,x0,kmax,N);
    zhist = Zgen(xnoisy,R,kmax,N);
    [Xkf,Pkf,ek] = KalmanFilterloc(Floc,r,sigmav2,sigmaw2,xhat0,P0loc,zhist,kmax,N,d,Aj,a,b,dT,xnoisy);
    time = 0:dT:kmax*dT;
    Xkftilde = Xkf-xnoisy;
    for i = 1:kmax
        MSEKF(i) = norm(Xkftilde(:,i));
    end
    
    %Plots
    if run == 1
    figure
    grid on
    l1 = plot(time,[x0(N+1:2*N) Xkf(N+1:2*N,:)],'r');
    grid on
    hold on
    grid on
    l2 = plot(time,[xhat0(N+1:2*N) xnoisy(N+1:2*N,:)],'b');
    grid on
    hold on
    l3 = plot(time,[0 rbar],'g');
    legend([l1(1), l2(1), l3(1)], {'Filtered Trajectories', 'Noisy Trajectories', 'Input'});
    title('Noisy Trajectories vs Filtered Trajectories')
    xlabel('time')
    ylabel('State/Input')
    axes('position',[.65 .175 .25 .25])
    box on
    plot(time(2000:2100),Xkf(N+1:2*N,2000:2100),'r',time(2000:2100),xnoisy(N+1:2*N,2000:2100),'b')
    axis tight
    
    figure
    plot(time,[0 MSEKF],'b');
    title('Mean-Square Error')
    xlabel('time')
    ylabel('MSE')
    end
    MSEKFAVG = mean(MSEKF);
    for i = 1:N
        e(i,run,:) = ek(i,1,:);
    end
end
alpha = .05;
r1 = chi2inv(alpha/2,Run)/Run;
r2 = chi2inv(1-alpha/2,Run)/Run;
sat = zeros(1,N);
for i = 1:N
    for j = 1:kmax
        ebar(i,j) = mean(e(i,:,j));
        if ebar(i,j)>r1 && ebar(i,j)<r2
            sat(i) = sat(i)+1;
        end
    end
    satperc(i) = sat(i)/kmax;
end
satperc
MSEKFAVG
eigenL = eig(L);
%neig = max(r(:,kmax))/(b*eigenL(2))