function [Oinf,converged]=max_pos_inv(Acl,X)

maxIterations=500;
Omega0 = X.minHRep(); % initialization
for i = 1:maxIterations
    % compute backward reachable set
    P = Pre_Aut(Acl,Omega0);
    % intersect with the state constraints
    P = P.intersect(Omega0).minHRep();
    if P == Omega0
        Oinf=Omega0;
        break
    else
        Omega0 = P;
    end
end
if i==maxIterations,
    converged=0;
else
    converged=1;
end
end



