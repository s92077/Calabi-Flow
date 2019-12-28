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

function G = GaussianCurvature(F, V, L, VB, VBno)
Fno = size(F,1);
Vno = size(V,1);
CosT = zeros(Fno,3);

% Compute dot product CosT = [ CosT_i, CosT_j, CosT_k ]
CosT(:,1) =  ( L(:,1).^2 + L(:,3).^2 - L(:,2).^2 ) ./ ( 2 * L(:,1) .* L(:,3) );
CosT(:,2) =  ( L(:,2).^2 + L(:,1).^2 - L(:,3).^2 ) ./ ( 2 * L(:,2) .* L(:,1) );
CosT(:,3) =  ( L(:,3).^2 + L(:,2).^2 - L(:,1).^2 ) ./ ( 2 * L(:,3) .* L(:,2) );

% Fix error
if ~isempty(find(CosT>1, 1)) || ~isempty(find(CosT<-1, 1))
    disp('Curvature Error!');
    disp(sum(CosT>1)+sum(CosT<-1));
    CosT(CosT>1)=1-eps;
    CosT(CosT<-1)=-1+eps;
end

% Compute W = [ theta_i, theta_j, theta_k ]
W(:,1) =  acos( CosT(:,1) );
W(:,2) =  acos( CosT(:,2) );
W(:,3) =  acos( CosT(:,3) );

% G is the weighted adjacency matrix
idx=(1:Fno)';
G = sum(sparse(F, [idx, idx, idx], W, Vno, Fno),2);
G = 2*pi - G;
if VBno > 0
    G(VB) = G(VB) - pi;
end
