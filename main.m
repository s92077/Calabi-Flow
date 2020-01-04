close all;
clear all;
[F,V,extra] = read_obj('hemisphere.obj');

delta=0.00001;

Vno=size(V,1);
[VB, VI, VBno] = BoundaryIndex(F);
E = sparse(F, F(:,[2,3,1]), true, Vno, Vno);
[I, R] = InvDistMetric(F, V, VB);
% [I, R] = ThurstonMetric(F, V, VB);
% [I, R] = TangentMetric(F, Vno);
% [I, R] = IsoscelesConstantAngleMetric(F, V, pi/2);
U=log(R);
iter=1;
while true
    L = CalculateLength(F, V, R, I);
    G = GaussianCurvature(F, V, L, VB, VBno);
%     G = GaussianCurvatureTangent(F, R, Vno, VB, VBno);
%     Gbar = zeros(size(G));
%     Gbar(VB,1) = 2*pi / VBno;
%     Gbar([45,46,1,82])=pi/2;
    Gbar = TargetCurvature(F, V, L, VB);
    fprintf('iter: %d %.9f\n',iter, max(abs(G-Gbar)));
    iter=iter+1;
    if max(abs(G-Gbar))<5e-8
        break;
    end
%     Ldual = DualLaplacian(F, V, R, I, L);
    EW = Hessian(F, Vno, R, L);
%     EW = RicciTangentOperator(F, V, R);
    EW(1,1)=EW(1,1)+1;
%     EW=EW+delta*speye(size(EW));
%     EW(sub2ind([Vno,Vno],1:Vno,1:Vno)) = EW(sub2ind([Vno,Vno],1:Vno,1:Vno))+1;
    U = U + EW\( Gbar - G );
%     U = U + delta * EW * ( Gbar - G );
    U = U - sum(U) / Vno;
    R = exp(U);
%     VL = Embedding(F, V, VB, L);
%     clf;
%     patch('Faces',F,'Vertices',VL,'FaceColor','green');
%     PlotMesh(F, VL);
%     patch(VL(VB,1),VL(VB,2),ones(VBno,1),'EdgeColor','interp','Marker','o','MarkerFaceColor','flat');
end
Z=euclidean_embed(F,U,I);
UV=[real(Z), imag(Z)];
PlotMesh(F,UV);

% figure
% subplot(1,2,1)
% PlotMesh(Ho.F, Ho.V, Ho.Vrgb);
% subplot(1,2,2)
% PlotMesh(F, UV, Ho.Vrgb);

AreaA=sum(abs(cross(reshape(UV(F,1),[size(F,1),3]),reshape(UV(F,2),[size(F,1),3]))),2);
AreaB=sum(cross(V(F(:,1),:)-V(F(:,2),:),V(F(:,1),:)-V(F(:,3),:)).^2,2);

L = CalculateLength(F, V, R, I);
CosT = zeros(size(F));
% Compute dot product CosT = [ CosT_i, CosT_j, CosT_k ]
CosT(:,1) =  ( L(:,1).^2 + L(:,3).^2 - L(:,2).^2 ) ./ ( 2 * L(:,1) .* L(:,3) );
CosT(:,2) =  ( L(:,2).^2 + L(:,1).^2 - L(:,3).^2 ) ./ ( 2 * L(:,2) .* L(:,1) );
CosT(:,3) =  ( L(:,3).^2 + L(:,2).^2 - L(:,1).^2 ) ./ ( 2 * L(:,3) .* L(:,2) );
AngleA=acos(CosT);
L = reshape(sqrt(sum((V(F(:,[2,3,1]),:) - V(F,:)).^2,2)), size(F));
CosT(:,1) =  ( L(:,1).^2 + L(:,3).^2 - L(:,2).^2 ) ./ ( 2 * L(:,1) .* L(:,3) );
CosT(:,2) =  ( L(:,2).^2 + L(:,1).^2 - L(:,3).^2 ) ./ ( 2 * L(:,2) .* L(:,1) );
CosT(:,3) =  ( L(:,3).^2 + L(:,2).^2 - L(:,1).^2 ) ./ ( 2 * L(:,3) .* L(:,2) );
AngleB=acos(CosT);
