function parameter = modelToGridUnits(model, parameter)
%MODELTOGRIDUNITS Convert the dimensionless parameters of the present model
%   (Averaging Theory notes, Section 2.3: unit periodic box [0,1]^2, length
%   scale L = box size, velocity scale V0, time scale L/V0, pressure scale
%   mu_f*V0/L) into the units the solvers work in.
%
%   The spectral grid (fourierGrid(1,N), fourierDiffOps) lives on [-1,1)^2,
%   i.e. box length Lgrid = 2 in grid units. One grid length unit is
%   therefore L/Lgrid, and every coefficient is rescaled by the power of L
%   it carries in eqs. (44)-(45):
%
%       kappa1 ~ L^2, F0 ~ L^2           ->  / Lgrid^2
%       S0 ~ L, Pc ~ L                    ->  / Lgrid
%       K ~ 1/L                           ->  * Lgrid
%       gamma1, eta_k (ratios to mu_f)    ->  unchanged
%       time t ~ L/V0                     ->  * Lgrid
%       pressure P ~ mu_f V0/L            ->  grid P = model P / Lgrid
%
%   Velocities are unaffected (both use V0).
%
%   Input  model     : struct with the model parameters (see mainModel.m)
%   Output parameter : the same struct the solvers receive, with the grid-
%                      unit coefficients added:
%       parameter.kappa1_0  drag coefficient at phi0   (kappa1 = kappa1_0*phi/phi0)
%       parameter.F0        active force               (F0*phi*d)
%       parameter.S0        active stress              (S0*div(phi*d d))
%       parameter.Pc        contact pressure           (grad(phi*Pc))
%       parameter.K0        Frank constant at phi0     (K = K0*phi/phi0)
%       parameter.gamma1    rotational friction        (constant, at phi0)
%       parameter.eta2      particle shear viscosity   (constant)
%       parameter.phi0      reference volume fraction

Lgrid = 2;   % fourierGrid(1,N): domain [-1,1)

parameter.Lgrid    = Lgrid;
parameter.phi0     = model.phi0;
parameter.kappa1_0 = model.kappa1_0 / Lgrid^2;
parameter.F0       = model.F0       / Lgrid^2;
parameter.S0       = model.S0       / Lgrid;
parameter.Pc       = model.Pc       / Lgrid;
parameter.K0       = model.K0       * Lgrid;
parameter.gamma1   = model.gamma1;
parameter.eta2     = model.eta2;
end
