%% compute localization length for Anderson model
% reference: Solid State Communications Vol. 67, No 3, pp 243
% reference: Z. Phys. B - Condensed Matter 62, 163-170 (1986)
N = 60000;
% nearest neighbour interaction
t = 1;
% disorder center
epsilon = 0.0;
% energy
E = 0.0;
% repeat times
NN = 1000;
lambda = zeros(NN,7);
% disorder width
for repeat =1:NN
    lambda_temp = zeros(1,7);
    indy = 1;
    for W = [1.0 1.5 2.0 2.5 3.0 3.5 4.0]
        indx = 2;
        % initial condition
        a = zeros(1,N);
        b = zeros(1,N);
        % random on-site potential
        V = (rand(1,N+1,1)-0.5)*W + epsilon;
        % initial condition
        a(1) = -t;
        b(1) = 0;
        
        temp = zeros(1,N);
        tol = 1;
        while abs(tol) >= 1e-3 && indx < N
            a(indx) = b(indx-1)-(V(indx)-E)*a(indx-1)/t;
            b(indx) = -a(indx-1);
            temp(indx) = indx/log(abs(a(indx)));
            tol = temp(indx);
            indx = indx + 1;
        end
        % get the index when a = 0
        M = indx - 2;
        % compute the localization length
        lambda_temp(indy) = mean(temp((M-500):M));
        indy = indy + 1;        
    end
    lambda(repeat,:) = lambda_temp;
end
loc = mean(lambda);
W = [1.0 1.5 2.0 2.5 3.0 3.5 4.0];
loglog(W,loc,'s')
% fitting the curve
log_W = log(W);
log_loc = log(loc);
P = polyfit(log_W,log_loc,1);
x = 1:0.01:4;
y = P(1)*log(x) + P(2);
hold on
% plot the fitting curve
loglog(x,exp(y),'k','linewidth',1.2)
% analytical results
%loc_analy = 105./W.^2;
%loglog(W,loc_analy,'*')
xlabel('disorder width W');
ylabel('\lambda')
%title('\lambda for Anderson Model, E = 0')
%legend('numerical','fitting')
