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

function Ldual = DualLaplacian(F, V, R, I, L)
Fno = size(F,1);
Vno = size(V,1);
Lratio = zeros(Fno,3);
LH = zeros(Fno,3);
CosT = zeros(Fno,3);

if size(V,2)==2
    V = [V, zeros(Vno,1)];
end

% Compute Rdiff = [Rk^2-Rj^2, Rj^2-Ri^2, Ri^2-Rk^2] 
Rdiff = [R(F(:,2)).^2 - R(F(:,1)).^2, R(F(:,3)).^2 - R(F(:,2)).^2, R(F(:,1)).^2 - R(F(:,3)).^2];

% Compute dot product Lratio = [ Lratio_ij, Lratio_jk, Lratio_ki ]
Lratio(:,1) =  (1 - Rdiff(:,1) ./ L(:,1).^2 ) / 2;
Lratio(:,2) =  (1 - Rdiff(:,2) ./ L(:,2).^2 ) / 2;
Lratio(:,3) =  (1 - Rdiff(:,3) ./ L(:,3).^3 ) / 2;

% Compute dot product LH = [ LHij, LHjk, LHki ]
LH(:,1) =  ( L(:,3).^2 + Rdiff(:,3) ) ./ ( 2 * L(:,1) .* L(:,3) );
LH(:,2) =  ( L(:,2).^2 + Rdiff(:,1) ) ./ ( 2 * L(:,2) .* L(:,1) );
LH(:,3) =  ( L(:,1).^2 + Rdiff(:,2) ) ./ ( 2 * L(:,3) .* L(:,2) );

% Compute dot product CosT = [ CosT_i, CosT_j, CosT_k ]
CosT(:,1) =  ( L(:,1).^2 + L(:,3).^2 - L(:,2).^2 ) ./ ( 2 * L(:,1) .* L(:,3) );
CosT(:,2) =  ( L(:,2).^2 + L(:,1).^2 - L(:,3).^2 ) ./ ( 2 * L(:,2) .* L(:,1) );
CosT(:,3) =  ( L(:,3).^2 + L(:,2).^2 - L(:,1).^2 ) ./ ( 2 * L(:,3) .* L(:,2) );

% W is the edge metric weight
W(:,1) = abs( (Lratio(:,1) .* CosT(:,1) - LH(:,1)) ./ sqrt( 1 - CosT(:,1).^2 ) );
W(:,2) = abs( (Lratio(:,2) .* CosT(:,2) - LH(:,2)) ./ sqrt( 1 - CosT(:,2).^2 ) );
W(:,3) = abs( (Lratio(:,3) .* CosT(:,3) - LH(:,3)) ./ sqrt( 1 - CosT(:,3).^2 ) );

% K is the weighted adjacency matrix
K = sparse(F, F(:,[2, 3, 1]), W, Vno, Vno);
K = K + K';

% L is the dual Laplacian operator
Ldual = K - diag( sum(K, 2) );
