%% compute localization length for 1D system
% disorder width
lambda = zeros(9,1);
ind = 1;
for W=1:9
    W
    % energy
    E = 0.0;
    N = 500000;
    % initial state
    %psi = [eye(M);zeros(M)];
    %R_all = zeros(2*M,2*M,N);
    T1 = eye(2);
    lambda1 = zeros(2,1);
    lambda2 = zeros(2,N);
    % number for renormalization
    nofororth = 10;
    
    for indx = 1:N
        for indy = 1:nofororth
            V = (rand-0.5)*W;
            T = [E-V -1;1 0];
            T1 = T*T1;
        end
        [Q,R] = qr(T1);
        D = diag(sign(diag(R)));
        Qunique = Q*D;
        Runique = D*R;
        T1 = Qunique;
        %R_all(:,:,ind) = Runique;
        for indz = 1:2
            lambda1(indz) = lambda1(indz) + log(Runique(indz,indz));
        end
        lambda2(:,indx) = lambda1/(nofororth*indx);
    end
    lambda1 = lambda1/(N*nofororth);
    lambda(ind) = lambda1(1);
    ind = ind + 1;
end
%figure
%plot(1:nofororth:N*nofororth,1./lambda2(1,:))
%ylim([100 500])
%xlabel('N')
%ylabel('$\lambda$','interpreter','latex')
%title('QR decomposition vs Transfer-Matrix')
%title('W = 1.0, E = 0.0')
W=1:9;
plot(W,lambda,'o')