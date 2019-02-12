%% exact diag one-dimensional Anderson localization
% problem
% system size
L = 5000;
lambda = zeros(9,1);
indx = 1;
% disorder width
for W = 1:9
    W
    temp = zeros(5,1);
    for ind = 1:5
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
%save('energy_W8.mat','W','eigval','L')