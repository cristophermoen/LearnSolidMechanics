% demo_q2.m
% Q2 — Full deformation field (u,v,w) and full 3x3 gradient matrix
% Uses uvw_prof_model(x,y,z,params); corner origin; y is vertical
clear; close all; clc;

%% ---- Parameters ----
params.Lx = 0.10;      % x-length (m)
params.H  = 0.10;      % y-thickness (vertical) (m)
params.Lz = 0.10;      % z-length (m)
inch = 0.0254;
params.Delta = 2.0 * inch;   % top compression at y=H (m)
params.A     = 0.02;        % bulging amplitude (m)

%% ---- 3D grid for visualization ----
Nx = 61; Nz = 61; Ny = 61;
xv = linspace(0, params.Lx, Nx);
zv = linspace(0, params.Lz, Nz);
yv = linspace(0, params.H , Ny);
[Xm, Zm, Ym] = meshgrid(xv, zv, yv);   % (x,z,y)

% Mid-planes for orthogonal slices
xs = params.Lx/2;  zs = params.Lz/2;  ys = params.H/2;

%% ---- Evaluate field & derivatives on grid ----
F = uvw_model(Xm, Ym, Zm, params);

%% ================= Figure A: u, v, w — 3D slices =================
figA = figure('Color','w','Position',[80 60 1300 430],'Renderer','opengl');
tA = tiledlayout(figA,1,3,'Padding','compact','TileSpacing','compact');

Ucells = {F.u, F.v, F.w};
Titles = {'u(x,y,z)', 'v(x,y,z)', 'w(x,y,z)'};

for k = 1:3
    ax = nexttile(tA);
    V = Ucells{k};
    slice(ax, Xm, Zm, Ym, V, xs, zs, ys);
    shading(ax,'interp'); axis(ax,'tight'); view(ax,3);
    xlabel(ax,'x (m)'); ylabel(ax,'z (m)'); zlabel(ax,'y (m)');
    title(ax, Titles{k}, 'Interpreter','tex'); colorbar(ax);

    % symmetric color scale
    cmax = max(abs(V(:)));
    if isempty(cmax) || ~isfinite(cmax) || cmax <= 0
        cmax = 1e-12;
    end
    caxis(ax, [-cmax cmax]);
end
sgtitle(tA,'Q2 — Deformation field (u, v, w) — 3D slices (y vertical)');

%% ============== Figure B: full 3×3 gradient matrix ==============
G = { F.ux, F.uy, F.uz; ...
      F.vx, F.vy, F.vz; ...
      F.wx, F.wy, F.wz };

Lbl = { '\partial u/\partial x','\partial u/\partial y','\partial u/\partial z'; ...
        '\partial v/\partial x','\partial v/\partial y','\partial v/\partial z'; ...
        '\partial w/\partial x','\partial w/\partial y','\partial w/\partial z' };

all_abs = cellfun(@(A) max(abs(A(:))), G);
cmaxB = max(all_abs(:));
if isempty(cmaxB) || ~isfinite(cmaxB) || cmaxB <= 0
    cmaxB = 1e-12;
end

figB = figure('Color','w','Position',[60 40 1300 1200],'Renderer','opengl');
tB = tiledlayout(figB,3,3,'Padding','compact','TileSpacing','compact');

for i = 1:3
    for j = 1:3
        ax = nexttile(tB);
        V = G{i,j};
        slice(ax, Xm, Zm, Ym, V, xs, zs, ys);
        shading(ax,'interp'); axis(ax,'tight'); view(ax,3);
        xlabel(ax,'x (m)'); ylabel(ax,'z (m)'); zlabel(ax,'y (m)');
        title(ax, Lbl{i,j}, 'Interpreter','tex'); colorbar(ax);
        caxis(ax, [-cmaxB cmaxB]);
    end
end
sgtitle(tB,'Q2 — Full gradient matrix \partial u_i/\partial x_j — 3D slices (y vertical)','Interpreter','tex');

%% ============== Part C: output deformation field at an arbitrary point ==============
x0 = params.Lx/3;
y0 = params.H/2;
z0 = params.Lz/4;

F0 = uvw_model(x0, y0, z0, params);

u0 = F0.u; v0 = F0.v; w0 = F0.w;

Grad = [F0.ux, F0.uy, F0.uz;
        F0.vx, F0.vy, F0.vz;
        F0.wx, F0.wy, F0.wz];

fprintf('Arbitrary point (x,y,z) = (%.4f, %.4f, %.4f) m\n', x0,y0,z0);
fprintf('Displacement: u = %.4e, v = %.4e, w = %.4e (m)\n', u0,v0,w0);
disp('Gradient of deformation field [du_i/dx_j]:');
disp(Grad);
