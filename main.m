close all;
clear all;
load('default.mat');
F=Ho.F;
V=Ho.V;
% PlotMesh(M.F, M.V);
% [F,V,extra] = read_obj('hemisphere.obj');

delta=0.001;

Vno=size(V,1);
[VB, VI, VBno] = BoundaryIndex(F);
E = sparse(F, F(:,[2,3,1]), true, Vno, Vno);
% [I, R] = InvDistMetric(F, V, VB);
[I, R] = ThurstonMetric(F, V, VB);
% [I, R] = TangentMetric(F, Vno);
U=log(R);
iter=1;
while true
    L = CalculateLength(F, V, R, I);
    G = GaussianCurvature(F, V, L, VB, VBno);
%     G = GaussianCurvatureTangent(F, R, Vno, VB, VBno);
    fprintf('iter: %d %f\n',iter, sum(abs(G)));
    iter=iter+1;
    if sum(abs(G))<0.1
        break;
    end
    Gbar = TargetCurvature(F, V, L, VB);
%     Gbar = zeros(size(G));
%     Gbar(VB,1) = 2*pi / VBno;
%     Gbar([45,46,1,82])=pi/2;
%     Ldual = DualLaplacian(F, V, R, I, L);
    EW = Hessian(F, Vno, R, L);
%     EW = RicciTangentOperator(F, V, R);
    EW(1,1)=EW(1,1)+1;
%     EW(sub2ind([Vno,Vno],1:Vno,1:Vno)) = EW(sub2ind([Vno,Vno],1:Vno,1:Vno))+1;
    U = U + EW\( Gbar - G );
%     U = U + Ldual * ( Gbar - G );
    U = U - sum(U) / Vno;
    R = exp(U);
%     VL = Embedding(F, V, VB, L);
%     clf;
%     patch('Faces',F,'Vertices',VL,'FaceColor','green');
%     PlotMesh(F, VL);
%     patch(VL(VB,1),VL(VB,2),ones(VBno,1),'EdgeColor','interp','Marker','o','MarkerFaceColor','flat');
end