%% exact diag one-dimensional Anderson localization with next-nearest-neighbor
% interaction, and compute localization length using Thouless formula
% problem
% system size
L = 3000;
% nearest neighbour
t1 = 1;
% next nearest neighbout
t2 = 0.01;
% energy
E = 0.0;
lambda = zeros(10,1);
errbar = zeros(10,1);
ind = 1;
% disorder center
epsilon = 0.0;

for W = linspace(1.5,2,10)
    A1 = zeros(5,1);
    for indx = 1:5
        % random potential
        Vp = (rand(L,1) - 0.5)*W + epsilon;
        V1 = t1*ones(L-1,1);
        V2 = t2*ones(L-2,1);
        
        H = diag(Vp)+diag(V1,-1)+diag(V1,1)+diag(V2,-2)+diag(V2,2);
        H(1,L) = t1;
        H(L,1) = t1;
        H(2,L) = t2;
        H(L,2) = t2;
        H(L-1,1) = t2;
        H(1,L-1) = t2;
        
        [v,d] = eig(H);
        eigval = diag(d);
        B = H(2:L,1:L-1);
        B = (-1)^(1+L)*det(B);
        A1(indx) = sum(log(abs(E-eigval)))/(L-1)-log(abs(B))/(L-1);
    end
    lambda(ind) = mean(A1);
    errbar(ind) = std(A1);
    ind = ind + 1
end
W = linspace(1.5,2,10);
W = W';
loc = 1./lambda;
loglog(W,loc,'ko')
log_W = log(W);
log_loc = log(loc);
P = polyfit(log_W,log_loc,1);
x = 1.5:0.01:2;
y = P(1)*log(x) + P(2);
hold on
loglog(x,exp(y),'k','linewidth',1.2)
xlabel('disorder width W');
ylabel('\lambda')
title('\lambda for NNN Model, t=1.0, T=1.0')
%legend('Thouless formula','y=-1.87x+4.5')
