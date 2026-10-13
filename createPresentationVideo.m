function createPresentationVideo(tt,CoarseMesh,kk,parameter,FrameStep) %#ok<INUSL> -- CoarseMesh kept for call-site compatibility, no longer used
% createPresentationVideo  Slimmed-down, large-font video for presentations.
% Shows only: volume fraction (phi), fluid velocity field, orientation field.
%
% Interpolation is FFT-based (fourierUpsample2D, O(N log N)) rather than
% the dense O(N^2) interpolation-matrix multiply -- same trigonometric
% interpolation kernel, see plotting_tool.m.
if nargin<5 || isempty(FrameStep)
    FrameStep = 5;
end

Ninterx=200;
Nintery=200;
Nx=parameter.Nx;
Ny=parameter.Ny;
[xx,~]=fourierGrid(1,Ninterx);
[yy,~]=fourierGrid(1,Nintery);
[XX,YY]=meshgrid(xx,yy);

idxPhi = 1:(Nx-2)*(Ny-2);
idxUf1 = (Nx-2)*(Ny-2)+1:2*Nx*Ny-2*Nx-2*Ny+4;
idxUf2 = 2*Nx*Ny-2*Nx-2*Ny+5:3*Nx*Ny-2*Nx-2*Ny+4;
idxTheta = 6*Nx*Ny-4*Nx-4*Ny+8:7*Nx*Ny-4*Nx-4*Ny+7;

% fixed color limits across the whole video (avoids flicker), taken
% directly from the raw coarse-grid values (cheap, no interpolation needed)
phiAll = kk(idxPhi,:);
phiMin = min(phiAll,[],'all');
phiMax = max(phiAll,[],'all');

speedAll = sqrt(kk(idxUf1,:).^2 + kk(idxUf2,:).^2);
speedMin = 0;
speedMax = max(speedAll,[],'all');

videoName = sprintf('PresentationVideo_%dx%d_S0_%g_K0_%g.mp4', Nx, Ny, ...
    parameter.model.S0, parameter.model.K0);
v = VideoWriter(videoName,'MPEG-4');
v.FrameRate = 10;
v.Quality = 100;
open(v)

fontSizeTitle = 26;
fontSizeAxes = 18;
fontSizeSup = 32;

fig = figure('units','normalized','outerposition',[0 0 1 1],'Color','w');

step = 8;

for i = 1:FrameStep:length(tt)

    t=tt(i);
    k=kk(:,i);
    psi=k(idxPhi);
    Uf1=k(idxUf1);
    Uf2=k(idxUf2);
    theta=k(idxTheta);

    figure(fig); clf;
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact')

    % --- Volume fraction ---
    nexttile
    phi_hires=reshape(fourierUpsample2D(psi, Ny-2, Nx-2, Nintery, Ninterx), Nintery, Ninterx);
    contourf(XX,YY,phi_hires,50,'LineColor','none')
    colormap(gca,turbo)
    clim([phiMin phiMax])
    cb = colorbar;
    cb.FontSize = fontSizeAxes;
    axis image
    xlim([-1 1])
    ylim([-1 1])
    title('Volume Fraction \phi','FontSize',fontSizeTitle)
    set(gca,'FontSize',fontSizeAxes)

    % --- Fluid velocity field ---
    nexttile
    U=reshape(fourierUpsample2D(Uf1, Ny, Nx, Nintery, Ninterx),Nintery,Ninterx);
    V=reshape(fourierUpsample2D(Uf2, Ny, Nx, Nintery, Ninterx),Nintery,Ninterx);
    contourf(XX,YY,sqrt(U.^2+V.^2),50,'LineColor','none')
    colormap(gca,turbo)
    clim([speedMin speedMax])
    cb = colorbar;
    cb.FontSize = fontSizeAxes;
    hold on
    quiver(XX(1:step:end,1:step:end), ...
           YY(1:step:end,1:step:end), ...
           U(1:step:end,1:step:end), ...
           V(1:step:end,1:step:end),'k','LineWidth',1.2,'AutoScaleFactor',1.5)
    axis image
    xlim([-1 1])
    ylim([-1 1])
    title('Fluid Velocity U_{fluid}','FontSize',fontSizeTitle)
    set(gca,'FontSize',fontSizeAxes)

    % --- Orientation field ---
    nexttile
    theta_field = reshape(fourierUpsample2D(theta, Ny, Nx, Nintery, Ninterx), Nintery, Ninterx);
    U = cos(theta_field);
    V = sin(theta_field);
    quiver(XX(1:step:end,1:step:end), ...
           YY(1:step:end,1:step:end), ...
           U(1:step:end,1:step:end), ...
           V(1:step:end,1:step:end),'k','LineWidth',1.2)
    axis image
    xlim([-1 1])
    ylim([-1 1])
    title('Orientation Field \theta','FontSize',fontSizeTitle)
    set(gca,'FontSize',fontSizeAxes)

    % t is in grid time units; model time unit is L/V0 = Lgrid grid units
    tModel = t/parameter.Lgrid;
    if isfield(parameter,'timeUnit_s') && ~isempty(parameter.timeUnit_s)
        sgtitle(sprintf('Time = %.2f s', tModel*parameter.timeUnit_s),'FontSize',fontSizeSup,'FontWeight','bold')
    else
        sgtitle(sprintf('Time = %.3f L/V_0', tModel),'FontSize',fontSizeSup,'FontWeight','bold')
    end

    drawnow

    frame = getframe(gcf);
    img = frame.cdata;

    [h, w, ~] = size(img);
    h = h - mod(h,2);
    w = w - mod(w,2);
    img = img(1:h, 1:w, :);

    writeVideo(v, img);
end
close(fig);
close(v)
end
