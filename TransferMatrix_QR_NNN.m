%% compute localization length for 1D NNN system
% disorder width
W = 1.0;
% energy
E = 0.0;
N = 500000;
t1 = 1.0;
t2 = 1.0;

% disorder center
epsilon = 0.0;
% initial state
%psi = [eye(M);zeros(M)];
%R_all = zeros(2*M,2*M,N);
T1 = eye(4);
lambda1 = zeros(4,1);
lambda2 = zeros(4,N);
% number for renormalization
nofororth = 10;

for indx = 1:N
    for indy = 1:nofororth
        V = (rand-0.5)*W + epsilon;
        T = [t1/t2 (E-V)/t2 t1/t2 1;1 0 0 0;0 1 0 0;0 0 1 0];
        T1 = T*T1;
    end
    [Q,R] = qr(T1);
    D = diag(sign(diag(R)));
    Qunique = Q*D;
    Runique = D*R;
    T1 = Qunique;
    %R_all(:,:,ind) = Runique;
    for indz = 1:4
        lambda1(indz) = lambda1(indz) + log(Runique(indz,indz));
    end
    lambda2(:,indx) = lambda1/(nofororth*indx);
end
lambda1 = lambda1/(N*nofororth);
%figure
plot(1:nofororth:N*nofororth,1./lambda2(2,:))
%ylim([100 500])
xlabel('N')
ylabel('$\lambda$','interpreter','latex')
%title('QR decomposition vs Transfer-Matrix')
title('W = 1.0, E = 0.1')