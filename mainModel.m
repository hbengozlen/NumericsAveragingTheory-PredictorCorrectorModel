clear all
clc
addpath('velocity', 'phi', 'theta', 'common')   % per-solver-family subdirectories
%% Predictor-corrector solver for the present two-phase model
% Dimensionless equations of the Averaging Theory notes, Section 2.3
% (eqs. 46-50), reduced as in Section 4.1: nu1 = nu2 = kappa2 = beta =
% Lambda = lambda = 0, eta_{1,3,4,5} = 0, so only the parameters below
% enter. Units: box side L, speed V0, time L/V0, pressure mu_f*V0/L.
% Constitutive closures: drag kappa1 = kappa1_0*phi/phi0, Frank constant
% K = K0*phi/phi0, all other coefficients constant.
%
% The solvers work on the [-1,1)^2 spectral grid (box length 2), so the
% model parameters are converted once by modelToGridUnits.m; all times
% below (T, dt, tt_model) are in model units L/V0.

%% Model parameters (dimensionless, notes Section 2.3)
% Default: Wolgemuth's Fig. 2 case mapped onto the model (notes Section
% 4.3, L = 400 um, V0 = 5 um/s => time unit L/V0 = 80 s). Reproduce with
% wolgemuthToModel.m, or set your own values directly here.
model.phi0     = 0.15;       % reference volume fraction phi_0
model.kappa1_0 = 9.796e4;    % drag at phi0            kappa1_0 = kappa1*L^2/mu_f      (45b, 81a)
model.F0       = 6.531e5;    % active force            F0 = F0*L^2/(mu_f*V0)           (44, 81a)
model.S0       = -6.531e3;   % active stress           S0 = S0*L/(mu_f*V0), <0 pushers (44, 81b)
model.Pc       = 3.265e3;    % contact pressure        Pc = Pc*L/(mu_f*V0)             (45c, 81b)
model.gamma1   = 0.8163;     % rotational friction     gamma1/mu_f                     (45a, 81b)
model.K0       = 0.4082;     % Frank constant at phi0  K/(mu_f*V0*L)                   (45b, 81c)
model.eta2     = 1;          % particle shear viscosity eta2/mu_f                      (45a)
                             %   (Wolgemuth used eta_b = 0 in his runs; kept 1 as in mainPC.m)

parameter.model      = model;      % keep the model-unit values for plotting/bookkeeping
parameter.timeUnit_s = 80;         % L/V0 in seconds, only used for video labels ([] = model units)
parameter            = modelToGridUnits(model, parameter);

parameter.PrefPt     = 100;   % pressure gauge mean (grid units)
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
dt_model         = 5e-5;     % = 1e-4 in mainPC.m's units (A/V0 = L/(2V0))
T_model          = 0.5;      % = T=1 in mainPC.m's units (40 s)
Nx               = 30;
Ny               = 30;
Theta0_amplitude = 1e-3;
Theta0_cutoff    = 10;

dt = dt_model * parameter.Lgrid;   % grid time units
T  = T_model  * parameter.Lgrid;

parameter.Nx = Nx;
parameter.Ny = Ny;

%% Grid and spectral operators
[parameter.x, ~] = fourierGrid(1, Nx);
[parameter.y, ~] = fourierGrid(1, Ny);
D = fourierDiffOps(Nx, Ny);
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
        t/parameter.Lgrid, mean(phiext_now), min(phiext_now), max(phiext_now), toc(stepWallTic));
    if min(phiext_now) < 0
        fprintf('phi_min=%.4e < 0, stopping simulation.\n', min(phiext_now));
        break;
    end
end

tt_model = tt / parameter.Lgrid;   % stored times in model units L/V0
% Note: kk holds velocities in units of V0 (same in both unit systems) and
% pressure in grid units -- multiply P by parameter.Lgrid for mu_f*V0/L.

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
