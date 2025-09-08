% demo_q5_prof.m
% Q5 – 3D strain field contours (Professor-based, y vertical)
clear; close all; clc;

% -------- Parameters --------
params.Lx = 0.10; params.H = 0.10; params.Lz = 0.10;
inch = 0.0254; params.Delta = 2.0*inch;
params.A = 0.002; 

% -------- 3D grid --------
Nx=41; Nz=41; Ny=41;
xv = linspace(0, params.Lx, Nx);
zv = linspace(0, params.Lz, Nz);
yv = linspace(0, params.H,  Ny);
[Xm, Zm, Ym] = meshgrid(xv, zv, yv);          % (x,z,y)

% -------- Strain tensor on grid --------
S = strain_tensor(Xm, Ym, Zm, params);

F3 = { ...
    S.eps_xx, S.eps_xy, S.eps_xz; ...
    S.eps_xy, S.eps_yy, S.eps_yz; ...
    S.eps_xz, S.eps_yz, S.eps_zz ...
    };

Labels = { ...
    '\epsilon_{xx}','\epsilon_{xy}','\epsilon_{xz}'; ...
    '\epsilon_{yx}','\epsilon_{yy}','\epsilon_{yz}'; ...
    '\epsilon_{zx}','\epsilon_{zy}','\epsilon_{zz}' ...
    };

% ---- draw strain cloud on 6 faces ----
% Field：'eq','exx','eyy','ezz','exy','exz','eyz'
plot_surface_strain(params, 'N', 61, 'Scale', 1.0, 'Field', 'eq');


% =====================================================================
% Figure 1: 3×3 Slices
% =====================================================================
xs = params.Lx/2; zs = params.Lz/2; ys = params.H/2;

fig1 = figure('Color','w','Position',[60 40 1400 1100],'Renderer','opengl');
t1 = tiledlayout(fig1,3,3,'Padding','compact','TileSpacing','compact');
for i=1:3
  for j=1:3
    V = F3{i,j};
    ax = nexttile(t1); 
    slice(ax, Xm, Zm, Ym, V, xs, zs, ys);
    shading(ax,'interp'); axis(ax,'tight'); view(ax,3);
    xlabel(ax,'x (m)'); ylabel(ax,'z (m)'); zlabel(ax,'y (m)');
    title(ax, ['Slice: ', Labels{i,j}], 'Interpreter','tex'); colorbar(ax);
    cmax = max(abs(V(:))); if isempty(cmax)||~isfinite(cmax)||cmax<=0, cmax=1e-12; end
    caxis(ax,[-cmax cmax]);
  end
end
sgtitle(t1,'Q5 – 3D Slices of Strain Components (y vertical)');