function model = wolgemuthToModel(w)
%WOLGEMUTHTOMODEL Map Wolgemuth's dimensional parameters onto the
%   dimensionless parameters of the model -- Table 1 of the paper
%   (Section 2.4). Only needed to reproduce/compare with Wolgemuth's case;
%   mainModel.m sets the model parameters directly.
%
%   Input struct w (SI units):
%     w.etaf      fluid viscosity mu_f                  [Pa s]
%     w.b         cell radius                           [m]
%     w.ell       cell length = force-dipole length     [m]
%     w.V0        swimming speed                        [m/s]
%     w.Dm        mixing diffusivity (contact pressure) [m^2/s]
%     w.KoverG1   K/gamma1 at phi0                      [m^2/s]  (1e-9 in the paper)
%     w.phi0      base volume fraction                  [-]
%     w.Lbox      side length of periodic box           [m]
%     w.pusher    true = pushers (S<0), false = pullers
%
%   Output: model struct with the fields mainModel.m expects
%   (shape/viscosity parameters set to the reference case: Lambda = 1,
%   lambda = 0, eta_k = nu = beta = kappa2 = 0).

L = w.Lbox; b = w.b;
sgn = 1 - 2*logical(w.pusher);                 % -1 pushers, +1 pullers

model.phi0      = w.phi0;
model.kappa1hat = 2*L^2/b^2;                   % kappa1 = 2 mu_f phi/b^2
model.F0        = 2*L^2/b^2;                   % F = 2 mu_f V0/b^2
model.S0        = sgn*2*w.ell*L/b^2;           % S = -2 mu_f l V0/b^2 (pushers)
model.Pc        = 2*w.Dm*L/(b^2*w.V0);         % Pc = 2 mu_f Dm/b^2
model.gamma1hat = w.ell^2/(6*b^2);             % gamma1 = mu_f l^2 phi/(6 b^2)
model.Khat      = w.KoverG1/(w.V0*L) * model.gamma1hat/w.phi0;   % K = Khat*phi^2
model.Lambda    = 1;  model.lambda = 0;  model.eta2 = 0;
model.eta1 = 0; model.eta3 = 0; model.eta4 = 0; model.eta5 = 0;
model.nu1  = 0; model.nu2  = 0; model.beta = 0; model.kappa2 = 0;
end
