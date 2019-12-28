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

function G = TargetCurvature(F, V, L, VB)
Fno = size(F,1);
Vno = size(V,1);
G = zeros(Vno,1);

% Compute the sparse edge length Ledge
Ledge = sparse(F,F(:,[2,3,1]),L,Vno, Vno);
Ledge = (Ledge + Ledge')/2;
Ledge(sub2ind([Vno, Vno], VB, circshift(VB,1))) = Ledge(sub2ind([Vno, Vno], VB, circshift(VB,1)))*2; 
Ledge(sub2ind([Vno, Vno], VB, circshift(VB,-1))) = Ledge(sub2ind([Vno, Vno], VB, circshift(VB,-1)))*2;

% Compute the boundary edge length Ledge
LB = Ledge(sub2ind([Vno, Vno], VB, circshift(VB,1)));

% Compute constant c = G_i / (l_(i-1 i) + l_(i-1 i)), sum(G_i)=2 pi
c  = pi / sum(LB);

% Compute target curvature G
G(VB) = c * (LB + circshift(LB,-1));