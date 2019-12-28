function W = TriangleVertices(U, V, L, Oref)
W = zeros([1,2]);
% Compute normal direction N
N = [V(2) - U(2), U(1) - V(1)];

% Compute dot product Cos = [ Lsum - Lratio_ij, Lsum - Lratio_jk, Lsum - Lratio_ki ]
Lsum = sum(L.^2) / 2;
Cratio = [Lsum - L(1).^2, Lsum - L(2).^2, Lsum - L(3).^2] / (L(1).^2);

% Compute W and Z
if L(3).^2 / L(1).^2 - Cratio(2).^2 < 0 
    disp(1);
end
W(1)=U(1) .* Cratio(2) + V(1) .* Cratio(3) + sign(sum((Oref - U) .* N, 2)) .* ...
    N(1) .* sqrt(L(3).^2 ./ L(1).^2 - Cratio(2).^2);
W(2)=U(2) .* Cratio(2) + V(2) .* Cratio(3) + sign(sum((Oref - U) .* N, 2)) .* ...
    N(2) .* sqrt(L(3).^2 ./ L(1).^2 - Cratio(2).^2);
end
