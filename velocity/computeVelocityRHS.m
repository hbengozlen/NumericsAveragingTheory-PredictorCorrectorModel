function rhs = computeVelocityRHS(phiext, d1, d2, F0, S0, Pc, boolactive, ...
    Ptilde, Np, RefPtIDXfine, D)
% RHS vectors (explicit terms) for solveVelocityFFT.m -- everything that
% does not multiply an unknown velocity/pressure, moved to the right:
%   fluid    (47):  -F0*phi*d + S0*div(phi*d d)       -> rhs  F0*phi*d - S0*div(phi*d d)
%   particle (48):  +F0*phi*d - grad(phi*Pc)          -> rhs -F0*phi*d + Pc*grad(phi)
% (lambda = 0, Pc constant). Ptilde carries the pressure gauge mean.
Ptilde_fine = D.liftP(Ptilde);
DxP0 = D.DXf(Ptilde_fine);
DyP0 = D.DYf(Ptilde_fine);

Fphi = F0*phiext;     % active force magnitude
Sphi = S0*phiext;     % active stress magnitude

if boolactive
    f1 =  Fphi.*d1 + (1-phiext).*DxP0 ...
       - D.DXf(Sphi.*d1.*d1) - D.DYf(Sphi.*d1.*d2);
    f2 =  Fphi.*d2 + (1-phiext).*DyP0 ...
       - D.DXf(Sphi.*d2.*d1) - D.DYf(Sphi.*d2.*d2);
    f3 = -Fphi.*d1 + phiext.*DxP0 + Pc*D.DXf(phiext);
    f4 = -Fphi.*d2 + phiext.*DyP0 + Pc*D.DYf(phiext);
else
    f1 = (1-phiext).*DxP0;
    f2 = (1-phiext).*DyP0;
    f3 = phiext.*DxP0 + Pc*D.DXf(phiext);
    f4 = phiext.*DyP0 + Pc*D.DYf(phiext);
end
f5 = zeros(Np,1);

f1(RefPtIDXfine) = 0;   % gauge pin RHS (matches velocityAv's pin)
f2(RefPtIDXfine) = 0;

rhs = [f1; f2; f3; f4; f5];
end
