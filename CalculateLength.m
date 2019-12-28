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

function L = CalculateLength(F, V, R, I)
Fno = size(F,1);
Vno = size(V,1);
L   = zeros(Fno,3);

% W is the length L = [ Lij, Ljk, Lki ];
L(:,1)=sqrt(R(F(:,1)).^2 + R(F(:,2)).^2 + 2 * R(F(:,1)) .* R(F(:,2)) ...
    .* I(sub2ind(size(I),F(:,1), F(:,2))));
L(:,2)=sqrt(R(F(:,2)).^2 + R(F(:,3)).^2 + 2 * R(F(:,2)) .* R(F(:,3)) ...
    .* I(sub2ind(size(I),F(:,2), F(:,3))));
L(:,3)=sqrt(R(F(:,3)).^2 + R(F(:,1)).^2 + 2 * R(F(:,3)) .* R(F(:,1)) ...
    .* I(sub2ind(size(I),F(:,3), F(:,1))));

