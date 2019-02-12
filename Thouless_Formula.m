load energy_W1.mat;
lambda = zeros(100,1);

for ind = 1:100
    lambda(ind) = sum(log(abs(eigval(:,ind))))/(5000-1);
end
%lambda = 1./lambda;
loc_length1 = mean(lambda);
err1 = std(lambda);


load energy_W2.mat;
lambda = zeros(100,1);

for ind = 1:100
    lambda(ind) = sum(log(abs(eigval(:,ind))))/(5000-1);
end
%lambda = 1./lambda;
loc_length2 = mean(lambda);
err2 = std(lambda);

load energy_W3.mat;
lambda = zeros(100,1);

for ind = 1:100
    lambda(ind) = sum(log(abs(eigval(:,ind))))/(5000-1);
end
%lambda = 1./lambda;
loc_length3 = mean(lambda);
err3 = std(lambda);


load energy_W4.mat;
lambda = zeros(100,1);

for ind = 1:100
    lambda(ind) = sum(log(abs(eigval(:,ind))))/(5000-1);
end
%lambda = 1./lambda;
loc_length4 = mean(lambda);
err4 = std(lambda);


load energy_W6.mat;
lambda = zeros(100,1);

for ind = 1:100
    lambda(ind) = sum(log(abs(eigval(:,ind))))/(5000-1);
end
%lambda = 1./lambda;
loc_length6 = mean(lambda);
err6 = std(lambda);

load energy_W7.mat;
lambda = zeros(100,1);

for ind = 1:100
    lambda(ind) = sum(log(abs(eigval(:,ind))))/(5000-1);
end
%lambda = 1./lambda;
loc_length7 = mean(lambda);
err7 = std(lambda);


load energy_W8.mat;
lambda = zeros(100,1);
for ind = 1:100
    lambda(ind) = sum(log(abs(eigval(:,ind))))/(5000-1);
end
%lambda = 1./lambda;
loc_length8 = mean(lambda);
err8 = std(lambda);

load energy_W9.mat;
lambda = zeros(100,1);
for ind = 1:100
    lambda(ind) = sum(log(abs(eigval(:,ind))))/(5000-1);
end
%lambda = 1./lambda;
loc_length9 = mean(lambda);
err9 = std(lambda);


loc_length = [loc_length1 loc_length2 loc_length3 loc_length4 loc_length6 loc_length7 loc_length8 loc_length9];
err = [err1 err2 err3 err4 err6 err7 err8 err9];

errorbar([1 2 3 4 6 7 8 9],loc_length,err,'sk','MarkerSize',1)
hold on
W = 1:0.01:10;
lambda = W.^2/105;
plot(W,lambda,'k-')
xlabel('disorder width W')
ylabel('1/\lambda')
legend('Thouless formula','W^2/105')
