function theta_new = solveThetaFFT(theta_old, phi, Uf1, Uf2, Up1, Up2, dt_step, parameter, D)
% Matrix-free, FFT-based twin of solveTheta.m -- same Backward Euler
% discretization, just applied without ever assembling solveTheta.m's
% dense Nfine x Nfine matrix. Solved with GMRES, RIGHT-preconditioned (see
% solveVelocityFFT.m for why right, not left: its relres is the TRUE
% residual, not a preconditioned stand-in) by the same operator with phi
% frozen to its mean (so its Laplacian term is exactly constant-
% coefficient) and Up1,Up2 frozen to their means.
%
% Warm-started from theta_old: under right preconditioning GMRES solves
% for y where x=thetaPinv(y), so passing theta_old as x0 directly would be
% wrong -- it needs mapping through thetaPinv's FORWARD (uninverted)
% operator first, so that thetaPinv(y0)~=theta_old. thetaPinv is exactly
% "divide by invDenom in Fourier space" (no lift/truncate -- theta already
% lives on the fine grid), so its exact forward map is "multiply by
% invDenom" instead, i.e. divide by thetaPinv's OWN divisor.
%
% D = fourierDiffOps(Nx,Ny), the same D the caller already builds once for
% solveVelocityFFT.m/solvePhiFFT.m -- Nx,Ny fixed for the whole run,
% passed in rather than recomputed here. phi is on the COARSE grid (same
% as solveTheta.m's CoarseMesh.interMatP*phi); D.liftP(phi) is its
% FFT-based equivalent (already the pattern used by solveVelocityFFT.m/
% computeVelocityRHS.m for the same lift).

Nx = parameter.Nx; Ny = parameter.Ny;
Nfine = Nx*Ny;
% Angular momentum balance (21d) with beta = nu2 = 0 gives
%   Omega3 = OmegaHat3_p + (1/gamma1) div(K grad theta)
% (+ Lambda*A_de, not yet included), inserted into the director equation
% (21e). The closures of Table 1 are evaluated here on the local phi:
%   K = Khat*phi^2,  gamma1 = gamma1hat*phi,
% and the elastic operator is applied in product-rule form
%   (1/gamma1) div(K grad v) = Kg.*lap(v) + gKx.*dv/dx + gKy.*dv/dy,
%   Kg = K/gamma1,  (gKx,gKy) = grad(K)/gamma1.
phiext = D.liftP(phi);
K      = parameter.Khat * phiext.^2;
gamma1 = parameter.gamma1hat * phiext;
Kg  = K ./ gamma1;
gKx = D.DXf(K) ./ gamma1;
gKy = D.DYf(K) ./ gamma1;

%% Frozen-mean preconditioner: K/gamma1 frozen to its mean (dropping the
% grad(K) terms) and Up1,Up2 frozen to their means makes the
% operator constant-coefficient, hence a per-mode scalar in Fourier space.
Kgbar = mean(Kg);
a = mean(Up1); b = mean(Up2);
invDenom = 1 ./ (1/dt_step + a*D.KX + b*D.KY - Kgbar*(D.KX2 + D.KY2));

ctx = struct('D', D, 'Up1', Up1, 'Up2', Up2, 'Kg', Kg, ...
    'gKx', gKx, 'gKy', gKy, ...
    'dt_step', dt_step, 'Nx', Nx, 'Ny', Ny, 'invDenom', invDenom);

rhs = computeThetaRHS(theta_old, Up1, Up2, dt_step, D);

y0 = real(ifft2(fft2(reshape(theta_old, Ny, Nx)) ./ invDenom));
y0 = y0(:);

tol = parameter.gmresTol;      % single shared value, set once in mainModel.m -- see its header comment
maxit = parameter.gmresMaxIt;  % single shared value, set once in mainModel.m -- see its header comment
AvPinv = @(y) thetaAv(thetaPinv(y, ctx), ctx);
[y, gmresFlag, gmresRelres, gmresIter] = gmres(AvPinv, rhs, min(parameter.gmresRestartScalar,Nfine), tol, maxit, [], [], y0);
theta_new = thetaPinv(y, ctx);

fprintf('  solveThetaFFT: gmres flag=%d relres=%.2e iters=%d/%d\n', ...
    gmresFlag, gmresRelres, gmresIter(1), gmresIter(2));
end
