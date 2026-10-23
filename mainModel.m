clear all
clc
addpath('velocity', 'phi', 'theta', 'common')   % per-solver-family subdirectories
%% Predictor-corrector solver for the present two-phase model
% Dimensionless equations (21) of "Linear stability of a two-phase model
% for moderately dense active suspensions" (Ben Gozlen, Wang, Oberlack).
% Units (Section 2.3): box side L, speed V0, time L/V0, pressure mu_f*V0/L.
%
% The code works directly in these units: the spectral grid is the unit
% periodic box (fourierDiffOps(Nx,Ny,1)), so all parameters, times (T, dt,
% tt), velocities and the pressure are exactly those of the paper. Grid
% coordinates parameter.x, parameter.y run over (0,1].

%% Model parameters: reference case, Table 1 of the paper
% Pusher B. subtilis following Wolgemuth (2008): cell radius b = 0.7 um,
% cell length = force-dipole length l = 4 um, V0 = 5 um/s,
% mu_f = 1e-3 Pa s, D_m = 1e-11 m^2/s, phi_bar = 0.15, L = 400 um
% (=> time unit L/V0 = 80 s). Closures: kappa1 = kappa1hat*phi,
% gamma1 = gamma1hat*phi, K = Khat*phi^2, everything else constant.
model.phi0      = 0.15;       % base volume fraction phi_bar
model.kappa1hat = 6.531e5;    % drag          kappa1 = kappa1hat*phi,  kappa1hat = 2L^2/b^2
model.F0        = 6.531e5;    % propulsion    F  = 2L^2/b^2
model.S0        = -6.531e3;   % active stress S  = -2 l L/b^2   (<0: pushers)
model.Pc        = 3.265e3;    % contact pressure Pc = 2 D_m L/(b^2 V0)
model.gamma1hat = 5.442;      % rot. friction gamma1 = gamma1hat*phi, gamma1hat = l^2/(6b^2)
model.Khat      = 18.14;      % Frank const.  K = Khat*phi^2, fixed by K/gamma1 = 1e-9 m^2/s
                              %   at phi_bar: Khat = (1e-9/(V0 L))*gamma1hat/phi_bar  (Er ~ 2.4)
model.Lambda    = 1;          % flow alignment, slender particles (Lambda = 1)
model.lambda    = 0;          % aspect ratio b/a -> 0 (slender limit)
model.eta2      = 0;          % particle viscosities eta_k = 0
% Zero in the reference case (Table 1) -- listed for completeness, the
% solvers do not contain these terms:
model.eta1 = 0; model.eta3 = 0; model.eta4 = 0; model.eta5 = 0;
model.nu1  = 0; model.nu2  = 0; model.beta = 0; model.kappa2 = 0;

parameter.model      = model;      % keep the model values for plotting/bookkeeping
parameter.timeUnit_s = 80;         % L/V0 in seconds, only used for video labels ([] = model units)
% Hand the model parameters to the solvers as they are -- the closures
% (kappa1hat*phi, gamma1hat*phi, Khat*phi^2) are evaluated inside the
% solvers on the local phi field.
f = fieldnames(model);
for i = 1:numel(f), parameter.(f{i}) = model.(f{i}); end
% NOTE: the solvers do not yet contain flow alignment (Lambda) or the
% elastic back-reaction on the particle momentum -- the parameters above
% are the Table 1 values, the equations still need those terms.
if parameter.Lambda ~= 0
    warning('mainModel:Lambda', 'Lambda = %g is set but not yet implemented in the solvers.', parameter.Lambda);
end

parameter.PrefPt     = 100;   % pressure gauge mean (units mu_f*V0/L)
parameter.boolactive = 1;

%% GMRES settings (see PredictorCorrector/mainPC.m for the reasoning behind each value)
parameter.useRightPrecond      = 1;     % right preconditioning (true residual)
parameter.gmresTol             = 1e-7;  % shared by all *FFT solvers
parameter.gmresMaxIt           = 100;   % shared outer-restart budget
parameter.gmresRestartVelocity = 100;   % Krylov size for solveVelocityFFT
parameter.gmresRestartScalar   = 30;    % Krylov size for phi/theta solvers

%% Run settings (model time units L/V0)
phi0             = model.phi0;
t0               = 0;
dt_model         = 5e-5;     % time step in L/V0
T_model          = 1;        % = 80 s
Nx               = 50;
Ny               = 50;
Theta0_amplitude = 1e-3;
Theta0_cutoff    = 10;

dt = dt_model;
T  = T_model;

parameter.Nx = Nx;
parameter.Ny = Ny;

%% Grid and spectral operators
parameter.Lbox = 1;                                       % unit box of the paper
[parameter.x, ~] = fourierGrid(parameter.Lbox/2, Nx);     % (-1/2, 1/2] ...
[parameter.y, ~] = fourierGrid(parameter.Lbox/2, Ny);
parameter.x = parameter.x + parameter.Lbox/2;             % ... shifted to (0, 1]
parameter.y = parameter.y + parameter.Lbox/2;
D = fourierDiffOps(Nx, Ny, parameter.Lbox);
CoarseMesh = [];   % only kept for plotting_tool/createPresentationVideo call compatibility

%% Initial conditions: uniform phi on the coarse grid, small random theta
Ncoarse = (Nx-2)*(Ny-2);
phi_n   = phi0 * ones(Ncoarse, 1);
theta_n = getTheta0(Nx, Ny, Theta0_cutoff, Theta0_amplitude);
theta_n = theta_n(:);
[Uf1_n, Uf2_n, Up1_n, Up2_n, P_n] = solveVelocityFFT(phi_n, theta_n, parameter, D);

% previous state for Step 4's extrapolated warm start (equal to n=0 at start)
Uf1_prev = Uf1_n; Uf2_prev = Uf2_n; Up1_prev = Up1_n; Up2_prev = Up2_n; P_prev = P_n;

%% Storage
t   = t0;
tt  = t0;
kk  = [phi_n; Uf1_n; Uf2_n; P_n; Up1_n; Up2_n; theta_n];

%% Profiling
perfCounter('solveVelocityFFT', 'reset');
timing = struct('solvePhi', 0, 'solveTheta', 0, 'solveVelocity', 0, ...
                 'solvePhiCN', 0, 'solveThetaCN', 0, 'nsteps', 0);
runTic = tic;

%% Time loop
while t < T
    stepWallTic = tic;
    dt_step = min(dt, T - t);

    % Step 1: predictor -- reuse current velocities
    Uf1_star = Uf1_n;  Uf2_star = Uf2_n;
    Up1_star = Up1_n;  Up2_star = Up2_n;

    % Step 2: phi half-step BE
    stepTic = tic;
    phi_half = solvePhiFFT(phi_n, Up1_star, Up2_star, dt_step/2, parameter, D);
    timing.solvePhi = timing.solvePhi + toc(stepTic);

    % Step 3: theta half-step BE
    stepTic = tic;
    theta_half = solveThetaFFT(theta_n, phi_n, Uf1_star, Uf2_star, Up1_star, Up2_star, dt_step/2, parameter, D);
    timing.solveTheta = timing.solveTheta + toc(stepTic);

    % Step 4: velocity solve with half-step state (extrapolated warm start)
    stepTic = tic;
    x0_h = 1.5*[Uf1_n; Uf2_n; Up1_n; Up2_n; P_n] ...
         - 0.5*[Uf1_prev; Uf2_prev; Up1_prev; Up2_prev; P_prev];
    [Uf1_h, Uf2_h, Up1_h, Up2_h, P_h] = solveVelocityFFT(phi_half, theta_half, parameter, D, x0_h);
    timing.solveVelocity = timing.solveVelocity + toc(stepTic);

    % Step 5: phi full CN step
    stepTic = tic;
    phi_n1 = solvePhiCNFFT(phi_n, Up1_h, Up2_h, dt_step, parameter, D);
    timing.solvePhiCN = timing.solvePhiCN + toc(stepTic);

    % Step 6: theta full CN step
    stepTic = tic;
    theta_n1 = solveThetaCNFFT(theta_n, phi_half, Uf1_h, Uf2_h, Up1_h, Up2_h, dt_step, parameter, D);
    timing.solveThetaCN = timing.solveThetaCN + toc(stepTic);

    % Step 7: velocity solve at full step (warm start from n and n+1/2)
    stepTic = tic;
    x0_n1 = 2*[Uf1_h; Uf2_h; Up1_h; Up2_h; P_h] - [Uf1_n; Uf2_n; Up1_n; Up2_n; P_n];
    [Uf1_n1, Uf2_n1, Up1_n1, Up2_n1, P_n1] = solveVelocityFFT(phi_n1, theta_n1, parameter, D, x0_n1);
    timing.solveVelocity = timing.solveVelocity + toc(stepTic);

    timing.nsteps = timing.nsteps + 1;

    Uf1_prev = Uf1_n; Uf2_prev = Uf2_n; Up1_prev = Up1_n; Up2_prev = Up2_n; P_prev = P_n;

    % Advance state
    phi_n   = phi_n1;
    theta_n = theta_n1;
    Uf1_n = Uf1_n1;  Uf2_n = Uf2_n1;
    Up1_n = Up1_n1;  Up2_n = Up2_n1;
    P_n   = P_n1;

    t   = t + dt_step;
    tt  = [tt, t];
    kk  = [kk, [phi_n; Uf1_n; Uf2_n; P_n; Up1_n; Up2_n; theta_n]];

    phiext_now = D.liftP(phi_n);
    fprintf('t=%.5f L/V0  phi_mean=%.4f  phi_min=%.4f  phi_max=%.4f  [step wall time: %.3f s]\n', ...
        t, mean(phiext_now), min(phiext_now), max(phiext_now), toc(stepWallTic));
    if min(phiext_now) < 0
        fprintf('phi_min=%.4e < 0, stopping simulation.\n', min(phiext_now));
        break;
    end
end

tt_model = tt;   % kept for compatibility: tt is already in units of L/V0
% kk: velocities in units of V0, pressure in units of mu_f*V0/L (paper units).

%% Timing report
totalTime = toc(runTic);
fprintf('\n=== Timing summary over %d steps (%.3f s total) ===\n', timing.nsteps, totalTime);
fprintf('%-14s %8.3f s  (%5.1f%%)\n', 'solvePhi',     timing.solvePhi,     100*timing.solvePhi/totalTime);
fprintf('%-14s %8.3f s  (%5.1f%%)\n', 'solveTheta',   timing.solveTheta,   100*timing.solveTheta/totalTime);
fprintf('%-14s %8.3f s  (%5.1f%%)\n', 'solveVelocity',timing.solveVelocity,100*timing.solveVelocity/totalTime);
fprintf('%-14s %8.3f s  (%5.1f%%)\n', 'solvePhiCN',   timing.solvePhiCN,   100*timing.solvePhiCN/totalTime);
fprintf('%-14s %8.3f s  (%5.1f%%)\n', 'solveThetaCN', timing.solveThetaCN, 100*timing.solveThetaCN/totalTime);

velocityFFTTiming = perfCounter('solveVelocityFFT', 'get', {'setup','Av','Pinv','calls','AvCalls','PinvCalls'});
fprintf('\n--- solveVelocityFFT breakdown (%d calls, %d Av calls, %d Pinv calls) ---\n', ...
    velocityFFTTiming.calls, velocityFFTTiming.AvCalls, velocityFFTTiming.PinvCalls);
gmresOverhead = timing.solveVelocity - velocityFFTTiming.setup - velocityFFTTiming.Av - velocityFFTTiming.Pinv;
fprintf('%-22s %8.3f s  (%5.1f%% of solveVelocity)\n', 'setup (pre-gmres)', velocityFFTTiming.setup, 100*velocityFFTTiming.setup/timing.solveVelocity);
fprintf('%-22s %8.3f s  (%5.1f%% of solveVelocity)\n', 'Av total',          velocityFFTTiming.Av,    100*velocityFFTTiming.Av/timing.solveVelocity);
fprintf('%-22s %8.3f s  (%5.1f%% of solveVelocity)\n', 'Pinv total',        velocityFFTTiming.Pinv,  100*velocityFFTTiming.Pinv/timing.solveVelocity);
fprintf('%-22s %8.3f s  (%5.1f%% of solveVelocity)\n', 'gmres overhead',    gmresOverhead,           100*gmresOverhead/timing.solveVelocity);
