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

function [I, Rmin] = TangentMetric(F, Vno)
% Set initial radius to be 1 and metric to be 0 to make the life simpler
Rmin = ones(Vno,1);
% W is the edge metric weight
W = ones(size(F));
% L is the initial Inversive Distance Circle Packing Metric
I = sparse(F, F(:,[2, 3, 1]), W, Vno, Vno);
I = (I + I') > 0;