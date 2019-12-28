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

function VL = Embedding(F, V, VB, L)
% C=cross(V(F(:,2),:)-V(F(:,1),:),V(F(:,3),:)-V(F(:,1),:));
% for i=1:size(F,1)
%     mArrow3(V(F(i,1),:),V(F(i,1),:)+C(i,:)./sqrt(sum(C(i,:).^2,2)));
% end
Fno = size(F,1);
Vno = size(V,1);
VL  = zeros(Vno, 2);

Lsparse = sparse(F, F(:,[2,3,1]), L, Vno, Vno);
Lsparse = (Lsparse + Lsparse') / 2;
Lsparse(sub2ind([Vno, Vno], VB, circshift(VB,1))) = Lsparse(sub2ind([Vno, Vno], VB, circshift(VB,1)))*2; 
Lsparse(sub2ind([Vno, Vno], circshift(VB,1), VB)) = Lsparse(sub2ind([Vno, Vno], VB, circshift(VB,-1)))*2;
FV      = sparse([1:Fno,1:Fno,1:Fno], F, ones(size(F)), Fno, Vno);

idx = [1];
Fembd = false(Fno,1);
Vembd = false(Vno,1);
c=1;


Vembd(F(1,1)) = true;
Vembd(F(1,2)) = true;
Vembd(F(1,3)) = true;
VL(F(1,1),:)=[ 0,0 ];
VL(F(1,2),:)=[ L(1,1),0 ];
VL(F(1,3),1)=( L(1,1)^2 + L(1,3)^2 - L(1,2)^2) / ( 2 * L(1,1) );
VL(F(1,3),2)=sqrt( L(1,3)^2 - VL(F(1,3),1)^2 );
patch('Faces',F(1,:),'Vertices',VL,'FaceColor','green');
while ~isempty(idx)
    t_idx=idx(1);
    t = F(t_idx, :);
    l = L(t_idx, :);
    for k=1:3
        t = circshift(t,1);
        l = circshift(l,1);
        v0=t(3);
        v1=t(1);
        v2=t(2);
        v3list=find( Lsparse(v1,:) & Lsparse(v2,:) );
        for v3=v3list
            tn = find( FV(:, v1) & FV(:, v2) & FV(:, v3));
            if Vembd(v3) || Fembd(tn)
                Fembd(tn)=true;
            else
                tn = tn(~Fembd(tn));
                idx = [idx, tn];

                find(F(tn,:),v1);
                find(F(tn,:),v2);
                VL(v3, :) = TriangleVertices(VL(v1,:), VL(v2,:), ...
                    [Lsparse(v1,v2), Lsparse(v2,v3), Lsparse(v3,v1)], VL(v0,:));
                Vembd(v3)=true;
%                 patch('Faces',F(tn,:),'Vertices',VL,'FaceColor','green');
            end
        end
    end
    if ~isempty(idx)
        idx=idx(2:end);
    end
    Fembd(t_idx)=true;
    c=c+1;
end
end

function W = TriangleVertices(U, V, L, Oref)
W = zeros([1,2]);
% Compute normal direction N
N = [V(2) - U(2), U(1) - V(1)];

% Compute dot product Cos = [ Lsum - Lratio_ij, Lsum - Lratio_jk, Lsum - Lratio_ki ]
Lsum = sum(L.^2) / 2;
Cratio = [Lsum - L(1).^2, Lsum - L(2).^2, Lsum - L(3).^2] / (L(1).^2);

% Compute W and Z
W(1)=U(1) .* Cratio(2) + V(1) .* Cratio(3) - sign(sum((Oref - U) .* N, 2)) .* ...
    N(1) .* sqrt(L(3).^2 ./ L(1).^2 - Cratio(2).^2);
W(2)=U(2) .* Cratio(2) + V(2) .* Cratio(3) - sign(sum((Oref - U) .* N, 2)) .* ...
    N(2) .* sqrt(L(3).^2 ./ L(1).^2 - Cratio(2).^2);
end
