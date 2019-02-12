%% this code is for Solid State Communications Vol. 67, No 3, pp 243
N = 5000000;
% disorder center
epsilon = 0.0;
% energy
E = 0.0;
% store the slope
p_slope = zeros(31,1);
indp = 1.0;
t = 1.0;
% nearest neighbout
for T = 1:0.01:1.3
    T
    indy = 1;
    lambda = zeros(10,1);
    for W = [0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0]
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
        while abs(tol) >= 1e-2 && indx <= N
            a(indx) = b(indx-1)-t*a(indx-1)/T;
            b(indx) = c(indx-1)-(V(indx+1)-E)*a(indx-1)/T;
            c(indx) = d(indx-1)-t*a(indx-1)/T;
            d(indx) = -a(indx-1);
            temp(indx) = indx/log(abs(a(indx)));
            tol = temp(indx);
            indx = indx + 1;
        end
        M = indx - 2;
        lambda(indy) = mean(temp((M-100):M));
        indy = indy + 1;        
    end
    W = [0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0];
    log_W = log(W);
    log_lambda = log(lambda);
    P = polyfit(log_W,log_lambda',1);
    p_slope(indp) = P(1);
    indp = indp + 1;
end