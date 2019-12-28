function EW = RicciTangentOperator(F, V, R)
    R = R(F);
    W = sqrt(R(:,1).*R(:,2).*R(:,3)./(R(:,1)+R(:,2)+R(:,3)));
    EW = sparse(F,F(:,[2,3,1]),repmat(W,[1,3])./( R + R(:,[2,3,1]) ) );
    EW = (EW + EW');
    EW = diag(sum(EW,2)) - EW;
end