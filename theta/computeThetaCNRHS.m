function rhs = computeThetaCNRHS(theta_old, Up1, Up2, dt_step, phiext, DXphiext, DYphiext, coeff, D)
% Explicit RHS for solveThetaCNFFT.m's Crank-Nicolson step:
%   rhs = theta_old/dt + 0.5*L_theta(theta_old) + f_src
% Literal transcription of solveThetaCN.m's dense
% rhs = (I/dt + 0.5*L_theta)*theta_old + f_src -- same L_theta = coeff*
% D_diff + D_adv (D_adv = -Up.grad(.)) as thetaAvCN.m, applied to the
% KNOWN theta_old instead of the unknown candidate v.
Dt1 = D.DXf(theta_old); Dt2 = D.DYf(theta_old);
D_diff_old = phiext.*(D.D2Xf(theta_old) + D.D2Yf(theta_old)) + DXphiext.*Dt1 + DYphiext.*Dt2;
D_adv_old  = -Up1.*Dt1 - Up2.*Dt2;
L_theta_old = coeff*D_diff_old + D_adv_old;

% Vorticity source OmegaHat3_p, sign as in eq. (50) of the notes (see
% computeThetaRHS.m); no zeta_f term in the present model.
f_src = 0.5*(D.DXf(Up2) - D.DYf(Up1));

rhs = theta_old/dt_step + 0.5*L_theta_old + f_src;
end
