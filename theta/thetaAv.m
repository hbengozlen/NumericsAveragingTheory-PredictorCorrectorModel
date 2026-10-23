function out = thetaAv(v, ctx)
% Matrix-free operator for solveThetaFFT.m: thetaAv(v) plays the role of
% A*v where A = I/dt + D_adv - D_diff, D_diff = (1/gamma1) div(K grad .) -- exact, not an
% approximation (same discretization as solveTheta.m's dense D_adv/D_diff,
% just applied without ever assembling a dense matrix). theta lives
% entirely on the fine grid already, so no lift/truncate round trip is
% needed here (unlike phiAv.m) -- just D.DXf/D.DYf/D.D2Xf/D.D2Yf directly.
D = ctx.D; Up1 = ctx.Up1; Up2 = ctx.Up2;
Kg = ctx.Kg; gKx = ctx.gKx; gKy = ctx.gKy; dt_step = ctx.dt_step;

Dv1 = D.DXf(v); Dv2 = D.DYf(v);
D_adv_v  = Up1.*Dv1 + Up2.*Dv2;
D_diff_v = Kg.*(D.D2Xf(v) + D.D2Yf(v)) + gKx.*Dv1 + gKy.*Dv2;   % (1/gamma1) div(K grad v)

out = v/dt_step + D_adv_v - D_diff_v;
end
