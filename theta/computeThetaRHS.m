function rhs = computeThetaRHS(theta_old, Up1, Up2, dt_step, D)
% Explicit RHS for solveThetaFFT.m: rhs = theta_old/dt + f_src, f_src the
% explicit vorticity source (no theta dependence):
%   f_src = OmegaHat3_p = 1/2*(dUp2/dx - dUp1/dy)
% -- sign as in the director equation (50) of the Averaging Theory notes
% with d = (cos theta, sin theta). Wolgemuth's fluid-vorticity term
% (zeta_f) is not part of the present model and has been removed.
f_src = 0.5*(D.DXf(Up2) - D.DYf(Up1));
rhs = theta_old/dt_step + f_src;
end
