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

function [I, Rmin] = ThurstonMetric(F,V,VB)
Fno = size(F,1);
Vno = size(V,1);
D = zeros(Fno,3);
Rv = zeros(Fno,3);

if size(V,2)==2
    V = [V, zeros(Vno,1)];
end

% Compute the cotangent weight
Vki = V(F(:,1),:) - V(F(:,3),:);
Vkj = V(F(:,2),:) - V(F(:,3),:);
Vij = V(F(:,2),:) - V(F(:,1),:);
Vc  = (V(F(:,1),:) + V(F(:,2),:) + V(F(:,3),:)) / 3;

% Compute D = [Dij, Djk, Dki] = [ Vij^2, Vjk^2, Vki^2 ]
D(:,1) =  sqrt( sum(Vij.^2, 2) );
D(:,2) =  sqrt( sum(Vkj.^2, 2) );
D(:,3) =  sqrt( sum(Vki.^2, 2) );

% Compute Rv = [Rv(i,jk), Rv(j,ki), Rv(k,ij)] 
Rv(:,1) =  (   D(:,1) - D(:,2) + D(:,3)  ) / 2;
Rv(:,2) =  (   D(:,1) + D(:,2) - D(:,3)  ) / 2;
Rv(:,3) =  ( - D(:,1) + D(:,2) + D(:,3)  ) / 2;

% Compute Rmin = [Rmin(i)] = [min_{jk} Rv(i,jk)]
idx   =  (1:Fno)';
R     =  sparse([idx, idx, idx], F, Rv, Fno, Vno);
Rmin  =  sum(R, 1)' ./ sum(R > 0, 1)';

% W is the edge metric weight including compensation for boundary edge
FVB = ismember(F,VB);
EB = FVB & FVB(:,[2,3,1]);
W(:,1)=( D(:,1).^2 - Rmin(F(:,1)).^2 - Rmin(F(:,2)).^2 ) ...
    ./ (2 * Rmin(F(:,1)) .* Rmin(F(:,2)));
W(:,2)=( D(:,2).^2 - Rmin(F(:,2)).^2 - Rmin(F(:,3)).^2 ) ...
    ./ (2 * Rmin(F(:,2)) .* Rmin(F(:,3)));
W(:,3)=( D(:,3).^2 - Rmin(F(:,3)).^2 - Rmin(F(:,1)).^2 ) ...
    ./ (2 * Rmin(F(:,3)) .* Rmin(F(:,1)));
W(W<0)=0;
W = W .* (1+EB);

% L is the initial Inversive Distance Circle Packing Metric
I = sparse(F, F(:,[2, 3, 1]), W, Vno, Vno);
I = (I + I') / 2;
% I(sub2ind([Vno, Vno], VB, circshift(VB,1))) = I(sub2ind([Vno, Vno], VB, circshift(VB,1)))*2; 
% I(sub2ind([Vno, Vno], VB, circshift(VB,-1))) = I(sub2ind([Vno, Vno], VB, circshift(VB,-1)))*2;
