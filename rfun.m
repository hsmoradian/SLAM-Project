function r = rfun(k,dT)
time = (k-1)*dT;
r = [3+0.5*sin(time);.5*cos(time);3*time/(1+time);time^2/(1+time^3);exp(-time/5)];