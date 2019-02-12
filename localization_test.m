%% this code is for Solid State Communications Vol. 67, No 3, pp 243
N = 20000000;
% nearest neighbout
t = 1;
% next nearest neighbour
T = 0.5;
% disorder center
epsilon = 0.0;
% energy
E = 0.0;
% disorder width
lambda = zeros(1,15);
indy = 1;
for W = [0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0]
    % initial condition
    a = zeros(1,N);
    b = zeros(1,N);
    c = zeros(1,N);
    d = zeros(1,N);
    V = (rand(1,N+1,1)-0.5)*W + epsilon;
    indx = 2;
    a(1) = t;
    b(1) = V(2)-E;
    c(1) = t;
    d(1) = T;
    temp = zeros(1,N);
    tol = 1;
    while abs(tol) >= 1e-2 && indx < N
        a(indx) = b(indx-1)-t*a(indx-1)/T;
        b(indx) = c(indx-1)-(V(indx+1)-E)*a(indx-1)/T;
        c(indx) = d(indx-1)-t*a(indx-1)/T;
        d(indx) = -a(indx-1);
        temp(indx) = indx/log(abs(a(indx)));
        tol = temp(indx);
        indx = indx + 1;
    end
    M = indx - 2
    lambda(indy) = mean(temp((M-100):M));
    indy = indy + 1;
    
end
W = [0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0];
%lambda = lambda/2;
loglog(W,lambda,'s')
log_W = log(W);
log_lambda = log(lambda);
P = polyfit(log_W,log_lambda,1);
x = 0.0:0.01:4;
y = P(1)*log(x) + P(2);
hold on
loglog(x,exp(y),'k','linewidth',1.2)
xlabel('disorder width W');
ylabel('\lambda')
title('\lambda for NNN Model, t=2, T=1.0')
%legend('numerical','y=-0.87x+2.74')
