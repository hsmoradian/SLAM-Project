clc
clear all

N = 5;
a = 1;
b = 10;
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
A = [On a*b*L;
     -I -a*I-b*L];
B = [On;a*I];
Sigmavm = [-a*b*Aj;b*Aj];
F = eye(2*N)+dT*A;
G = dT*B;
Gamma = dT*Sigmavm;
sigmav2 = .0001;
sigmaw2 = .1;
v = sigmav2*ones(1,N);
w = sigmaw2*ones(1,N);
Q = diag(v);
R = diag(w);
H = [On I];
intervalx0 = 10;
x0 = [zeros(N,1);intervalx0*rand(N,1)];
xhat0 = [zeros(N,1);intervalx0*rand(N,1)];
P0 = 10*eye(2*N);
tmax = 100;
kmax = tmax/dT;
for i = 1:kmax
    r(:,i) = rfun(i,dT);
    rbar(1,i) = mean(r(:,i));
end
Run = 1;
for run = 1:Run
    xnoisy = Xgennoisy(F,G,Gamma,r,Q,x0,kmax,N);
    zhist = Zgen(xnoisy,R,kmax,N);
    [Xkf,Pbarkf,Pkf] = KalmanFilter(F,G,r,Gamma,H,Q,R,xhat0,P0,zhist,kmax);
    time = 0:dT:kmax*dT;
    Xkftilde = Xkf-xnoisy;
    for i = 1:kmax
        MSEKF(i) = norm(Xkftilde(:,i));
    end
    %Plots
    
    %{
    %}
    
    figure
    plot(time,[xhat0 xnoisy],'b')
    title('Noisy')
    hold on
    plot(time,[x0 Xkf],'r')
    title('Kalman')
    figure
    plot(time,[0 MSEKF])
    title('Kalman Error')
    MSEKFAVG = mean(MSEKF);
    
    e(run,:) = egen(Xkftilde,Pkf,kmax);
end
alpha = .01;
r1 = chi2inv(alpha/2,2*Run*N)/Run;
r2 = chi2inv(1-alpha/2,2*Run*N)/Run;
sat = 0;
for i = 1:kmax
    ebar(i) = mean(e(:,i));
    if ebar(i)>r1 && ebar(i)<r2
        sat = sat+1;
    end
end
satperc = sat/kmax
