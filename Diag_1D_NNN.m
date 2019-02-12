tic
%% exact diag one-dimensional Anderson localization
% problem
% system size
L = 2000;
% nearest neighbour
t1 = 1.0;
% next nearest neighbout
t2 = 0.01;
% energy
E = 0.0;
% repeat number
M = 10;
lambda = zeros(11,1);
errbar = zeros(11,1);
ind = 1;
% disorder center
epsilon = 0.0;
%for ind = 1:10
for W = [0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0]
    A1 = zeros(M,1);
    for indx = 1:M
        % random potential
        Vp = (rand(L,1) - 0.5)*W + epsilon;
        V1 = t1*ones(L-1,1);
        V2 = t2*ones(L-2,1);
        
        H = diag(Vp) + diag(V1,-1) + diag(V1,1) + diag(V2,-2) + diag(V2,2);
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

W = [0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0]';
loc = 1./lambda;
%loc_length = W.^2/105;
%loc_length = 105./W.^2;
loglog(W,loc,'k*')
log_W = log(W);
log_loc = log(loc);
P = polyfit(log_W,log_loc,1);
x = 0.5:0.01:1.0;
y = P(1)*log(x) + P(2);
hold on
loglog(x,exp(y),'k','linewidth',1.2)
xlabel('disorder width W');
ylabel('\lambda')
title('\lambda for NNN Model, t=1.0,T=0.2')
%legend('Thouless formula','y=-1.87x+4.5')
%plot(W,loc_length,'-k')
%hold on
%errorbar(1:10,lambda,errbar,'sk','MarkerSize',1)
%plot(1:10,1./lambda,'o')
%legend('W^2/105','Thouless formula')
%save('t1.mat')
toc