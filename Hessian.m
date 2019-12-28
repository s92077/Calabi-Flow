%% LaplaceBeltrami
%  Compute the discrete Laplace-Beltrami operator.
%  Please refer to [1] for more details.
%  [1] M.-H. Yueh, W.-W. Lin, C.-T. Wu, and S.-T. Yau, 
%      An efficient energy minimization for conformal parameterizations, 
%      J Sci Comput (2017). doi:10.1007/s10915-017-0414-y
%
%% Syntax
%   L = LaplaceBeltrami(F, V)
%
%% Description
%  F  : double array, nf x 3, faces of mesh
%  V  : double array, nv x 3, vertices of mesh
% 
%  L  : double array, nv x nv, matrix of Laplaci-Beltrami operator of mesh
%
%% Contribution
%  Author : Mei-Heng Yueh
%  Created: 2016/09/06
% 
%  Copyright 2016 Mei-Heng Yueh
%  http://scholar.harvard.edu/yueh

function W = Hessian(F, Vno, R, L)
Fno = size(F,1);
C   = zeros(Fno,3);

T=@(li,rj,rk) 0.5*( li.*li+rj.*rj-rk.*rk );

K = zeros(Fno,9);
W = sparse(Vno, Vno);

Lsum = sum(L,2) / 2;
A    = sqrt(Lsum .* (Lsum - L(:,1)) .* (Lsum - L(:,2)) .* (Lsum - L(:,3)));
% Compute dot product CosT = [ CosT_i, CosT_j, CosT_k ]
C(:,1) =  ( L(:,1).^2 + L(:,3).^2 - L(:,2).^2 ) ./ ( 2 * L(:,1) .* L(:,3) );
C(:,2) =  ( L(:,2).^2 + L(:,1).^2 - L(:,3).^2 ) ./ ( 2 * L(:,2) .* L(:,1) );
C(:,3) =  ( L(:,3).^2 + L(:,2).^2 - L(:,1).^2 ) ./ ( 2 * L(:,3) .* L(:,2) );

for t=1:Fno
    Lvec = L(t,[2,3,1]);
    Theta = [-1      C(t,3)  C(t,2) ; ...
             C(t,3)  -1      C(t,1) ; ...
             C(t,2)  C(t,1)  -1     ];
    D     = [0                              T(L(t,2),R(F(t,2)),R(F(t,3)))  T(L(t,2),R(F(t,3)),R(F(t,2))) ; ...
             T(L(t,3),R(F(t,1)),R(F(t,3)))  0                              T(L(t,3),R(F(t,3)),R(F(t,1))) ; ...
             T(L(t,1),R(F(t,1)),R(F(t,2)))  T(L(t,1),R(F(t,2)),R(F(t,1)))  0                      ];
    K(t,:) = reshape(diag(Lvec) * Theta * diag(1./Lvec) * D / A(t), [1,9]);
    B = diag(Lvec) * Theta * diag(1./Lvec) * D / A(t);
%     if sum(B(1,:),2)~=0 || sum(B(2,:),2)~=0 || sum(B(3,:),2)~=0
%         disp(1);
%     end
end
W = sparse([F F F],[repmat(F(:,1),[1,3]) repmat(F(:,2),[1,3]) repmat(F(:,3),[1,3])],K, Vno, Vno)/2;