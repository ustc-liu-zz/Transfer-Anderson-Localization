%% compute localization length for 1D NNN system
% disorder width
%W = 1.0;
% energy
E = 0.0;
N = 500000;
t1 = 0.1;
t2 = 1.0;

% disorder center
epsilon = 0.0;
% initial state
%psi = [eye(M);zeros(M)];
%R_all = zeros(2*M,2*M,N);
lambda = zeros(1,1);
ind = 1;
for W = 1:0.5:4.0
    T1 = eye(4);
    temp1 = zeros(4,1);
    temp2 = zeros(4,N);
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
            temp1(indz) = temp1(indz) + log(Runique(indz,indz));
        end
        temp2(:,indx) = temp1/(nofororth*indx);
    end
    temp1 = temp1/(N*nofororth);
    lambda(ind) = 1/temp1(2);
    ind = ind + 1;
end
%figure
%plot(1:nofororth:N*nofororth,1./lambda(2,:))
%ylim([100 500])
%xlabel('N')
%ylabel('$\lambda$','interpreter','latex')
%title('QR decomposition vs Transfer-Matrix')
%title('W = 1.0, E = 0.1')