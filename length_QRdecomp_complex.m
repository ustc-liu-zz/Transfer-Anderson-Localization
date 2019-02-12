%% test the difference bewteen QR factoration and transfer-matrix method in two-dimensional system
% disorder width
W1 = sqrt(3);
W2 = 0;
% energy
E = 0.5;
% system width
M = 16;
N = 100000;
% initial state
%psi = [eye(M);zeros(M)];
%R_all = zeros(2*M,2*M,N);
T1 = eye(2*M);
lambda1 = zeros(2*M,1);
lambda2 = zeros(2*M,N);
% number for renormalization
nofororth = 10;

for indx = 1:N
    for indy = 1:nofororth
        V = (rand(M,1)-0.5)*W1 + (rand(M,1)-0.5)*W2*1i;
        H = diag(V)+diag(ones(M-1,1),1)+diag(ones(M-1,1),-1);
        H(1,M) = 1;
        H(M,1) = 1;
        T = [-E*eye(M)+H eye(M);-eye(M) zeros(M)];
        T1 = T*T1;
    end
    [Q,R] = qr(T1);
    D = diag(sign(diag(R)));
    Qunique = Q*D;
    Runique = D*R;
    T1 = Qunique;
    %T1 = R;
    %R_all(:,:,ind) = Runique;
    for indz = 1:2*M
        lambda1(indz) = lambda1(indz) + log(abs(R(indz,indz)));
    end
    lambda2(:,indx) = lambda1/(nofororth*indx);
end
lambda1 = lambda1/(N*nofororth);
%figure
plot(1:nofororth:N*nofororth,1./lambda2(16,:),'-.','linewidth',1.5)
ylim([-500 800])
xlabel('N')
ylabel('$\lambda$','interpreter','latex')
title('QR decomposition')