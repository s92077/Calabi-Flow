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

function L = LaplaceBeltrami(F,V)
Fno = size(F,1);
Vno = size(V,1);
W = zeros(Fno,3);

if size(V,2)==2
    V = [V, zeros(Vno,1)];
end

% Compute the cotangent weight
Vki = V(F(:,1),:) - V(F(:,3),:);
Vkj = V(F(:,2),:) - V(F(:,3),:);
Vij = V(F(:,2),:) - V(F(:,1),:);

% Compute W = [Wij, Wjk, Wki] = [ cot(theta_k), cot(theta_i), cot(theta_j)]
W(:,1) =  0.5*sum(Vki.*Vkj, 2) ./ sqrt( sum(cross(Vki, Vkj).^2, 2) );
W(:,2) = -0.5*sum(Vij.*Vki, 2) ./ sqrt( sum(cross(Vij, Vki).^2, 2) );
W(:,3) =  0.5*sum(Vkj.*Vij, 2) ./ sqrt( sum(cross(Vkj, Vij).^2, 2) );

% K is the weighted adjacency matrix
K = sparse(F, F(:,[2, 3, 1]), W, Vno, Vno);
K = K + K';

% L is the discrete Laplaci-Beltrami operator
L = diag( sum(K, 2) ) - K;
