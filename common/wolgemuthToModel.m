function model = wolgemuthToModel(w)
%WOLGEMUTHTOMODEL Map Wolgemuth's dimensional parameters onto the
%   dimensionless parameters of the present model -- Table 1 and eq. (81)
%   of the Averaging Theory notes (Section 4.2-4.3). Only needed to
%   reproduce/compare with Wolgemuth's case; mainModel.m sets the model
%   parameters directly.
%
%   Input struct w (SI units):
%     w.etaf      fluid viscosity                 [Pa s]
%     w.etab      particle viscosity eta_b        [Pa s]
%     w.R         bacterial radius                [m]
%     w.Lb        bacterial length                [m]
%     w.b         active-stress dipole length     [m]
%     w.V0        swimming speed                  [m/s]
%     w.Db        diffusivity in contact pressure [m^2/s]
%     w.KoverZr   K/zeta_r as in eq. (81c)         [m^2/s]  (1e-9 in the notes;
%                 mainPC.m's elasPara corresponds to elasPara*6*R^2*phi0/(etaf*Lb^2))
%     w.phi0      reference volume fraction       [-]
%     w.Lbox      side length of periodic box     [m]
%     w.pusher    true = pushers (S0<0), false = pullers
%
%   Output: model struct with the fields mainModel.m expects.

L = w.Lbox;
sgn = 1 - 2*logical(w.pusher);           % -1 pushers, +1 pullers

model.phi0     = w.phi0;
model.kappa1_0 = 2*w.phi0*L^2 / w.R^2;                     % (81a)
model.F0       = 2*L^2 / w.R^2;                            % (81a)
model.S0       = sgn * 2*w.b*L / w.R^2;                    % (81b)
model.Pc       = 2*w.Db*L / (w.R^2*w.V0);                  % (81b)
model.gamma1   = w.Lb^2*w.phi0 / (6*w.R^2);                % (81b)
model.K0       = w.KoverZr * model.gamma1 / (w.V0*L);      % (81c)
model.eta2     = w.etab / w.etaf;                          % Table 1: eta_2 = eta_b
end
