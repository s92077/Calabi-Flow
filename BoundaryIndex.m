%% BoundaryIndex
%  Compute the boundary index of a simply connected open mesh.
%
%% Syntax
%   [VB, VI, VBno] = BoundaryIndex(F)
%
%% Description
%  F    : double array, nf x 3, faces of mesh
%
%  VB   : double array, VBno x 1, index of boundary vertices
%  VI   : double array, VIno x 1, index of interior vertices
%  VBno : double, number of boundary vertices
%
%% Contribution
%  Author : Mei-Heng Yueh
%  Created: 2016/09/06
% 
%  Copyright 2016 Mei-Heng Yueh
%  http://scholar.harvard.edu/yueh

function [VB, VI, VBno] = BoundaryIndex(F)
Vno = max(max(F));
Gvv = sparse(F, F(:,[2 3 1]), 1, Vno, Vno);
Gb  = Gvv - Gvv.';
[Bi, Bj] = find( Gb == 1 );
VBno = numel(Bi);
VB   = zeros(VBno, 1);
[bd_vertex, bd_ind] = min(Bi);
for jj = 1:VBno
    VB(jj) = bd_vertex;
    bd_ind = find( Bi == Bj(bd_ind) );
    bd_vertex = Bi(bd_ind);
    if ( bd_vertex == VB(1) ) && ( jj ~= VBno )
        VBno = jj;
        VB(jj+1:end) = [];
        break
    end
end
if nargout > 1
    VI = setdiff((1:Vno).', VB);
end