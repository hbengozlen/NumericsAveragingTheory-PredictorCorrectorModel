function theta_new = solveThetaCNFFT(theta_old, phi, Uf1, Uf2, Up1, Up2, dt_step, parameter, D)
% Matrix-free, FFT-based twin of solveThetaCN.m -- same Crank-Nicolson
% discretization (including its D_adv sign convention -- see
% solveThetaCN.m's comment), just applied without ever assembling
% solveThetaCN.m's dense Nfine x Nfine matrix. Solved with GMRES, RIGHT-
% preconditioned (see solveVelocityFFT.m for why right, not left) by the
% same operator with phi frozen to its mean and Up1,Up2 frozen to their
% means.
%
% Warm-started from theta_old via thetaPinvCN's forward (uninverted)
% operator -- see solveThetaFFT.m's comment for why (same reasoning, CN's
% invDenomCN in place of BE's invDenom).
%
% D = fourierDiffOps(Nx,Ny), the same D the caller already builds once for
% solveVelocityFFT.m/solvePhiFFT.m. phi is on the COARSE grid;
% D.liftP(phi) is the FFT-based equivalent of solveThetaCN.m's
% CoarseMesh.interMatP*phi.

Nx = parameter.Nx; Ny = parameter.Ny;
Nfine = Nx*Ny;
% Angular momentum (49) with beta=Lambda=nu2=0 gives
%   Omega3 = OmegaHat3_p + (1/gamma1) div(K grad theta),
% inserted into the director equation (50). K = K0*phi/phi0 (linear in
% phi), gamma1 constant -- this reproduces the previous Wolgemuth-based
% operator exactly, now written with the model's K0, gamma1 (grid units,
% modelToGridUnits.m).
coeff = parameter.K0 / (parameter.gamma1 * parameter.phi0);   % div(K grad)/gamma1 = coeff*div(phi grad)

phiext = D.liftP(phi);
DXphiext = D.DXf(phiext);
DYphiext = D.DYf(phiext);

%% Frozen-mean preconditioner: same reasoning as solveThetaFFT.m's, just
% with the (I/dt - 0.5*L_theta) form Crank-Nicolson puts on the implicit
% operator instead of BE's (I/dt + D_adv - coeff*D_diff).
phibar = mean(phiext);
a = mean(Up1); b = mean(Up2);
invDenomCN = 1 ./ (1/dt_step - 0.5*(coeff*phibar*(D.KX2 + D.KY2) - a*D.KX - b*D.KY));

ctx = struct('D', D, 'Up1', Up1, 'Up2', Up2, 'phiext', phiext, ...
    'DXphiext', DXphiext, 'DYphiext', DYphiext, 'coeff', coeff, ...
    'dt_step', dt_step, 'Nx', Nx, 'Ny', Ny, 'invDenomCN', invDenomCN);

rhs = computeThetaCNRHS(theta_old, Up1, Up2, dt_step, ...
    phiext, DXphiext, DYphiext, coeff, D);

y0 = real(ifft2(fft2(reshape(theta_old, Ny, Nx)) ./ invDenomCN));
y0 = y0(:);

tol = parameter.gmresTol;      % single shared value, set once in mainModel.m -- see its header comment
maxit = parameter.gmresMaxIt;  % single shared value, set once in mainModel.m -- see its header comment
AvPinv = @(y) thetaAvCN(thetaPinvCN(y, ctx), ctx);
[y, gmresFlag, gmresRelres, gmresIter] = gmres(AvPinv, rhs, min(parameter.gmresRestartScalar,Nfine), tol, maxit, [], [], y0);
theta_new = thetaPinvCN(y, ctx);

fprintf('  solveThetaCNFFT: gmres flag=%d relres=%.2e iters=%d/%d\n', ...
    gmresFlag, gmresRelres, gmresIter(1), gmresIter(2));
end
