%% exact diag one-dimensional Anderson localization and compute localization
% length using Thouless formula
% problem
% system size
L = 5000;
lambda = zeros(7,1);
indx = 1;
% disorder width
for W = 1:0.5:4
    W
    temp = zeros(10,1);
    for ind = 1:10
        % random potential
        Vp = (rand(L,1) - 0.5)*W;
        
        H = diag(Vp) + diag(ones(L-1,1),-1)+diag(ones(L-1,1),1);
        H(1,L) = 1;
        H(L,1) = 1;
        
        [v,d] = eig(H);
        eigval = diag(d);
        temp(ind) = sum(log(abs(eigval)))/(5000-1);
        
    end
    lambda(indx) = mean(1./temp);
    indx = indx + 1;
end
W = [1.0 1.5 2.0 2.5 3.0 3.5 4.0];
loglog(W,lambda,'*')
% fitting the curve
log_W = log(W);
log_loc = log(lambda);
P = polyfit(log_W,log_loc',1);
x = 1:0.01:4;
y = P(1)*log(x) + P(2);
hold on
% plot the fitting curve
loglog(x,exp(y),'-k','linewidth',1.2)
%save('energy_W8.mat','W','eigval','L')