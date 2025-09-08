clear; close all; clc;

% ========= Parameters =========
params.Lx = 0.10;    % m
params.Ly = 0.10;    % m
params.H  = 0.10;    % m
inch = 0.0254;
params.Delta = 2.0*inch;   % Top surface pressed down 2 inches
params.A  = 0.02;         % Vertical flexure amplitude
params.B  = 0.01;         % Lateral displacement amplitude

% ========= 1. Strain tensor at any point =========
x0 = 0.00; y0 = 0.00; z0 = params.H/2;  % Example: Center of medium thickness
S0 = strain_tensor(x0,y0,z0,params);

E = [ S0.eps_xx, S0.eps_xy, S0.eps_xz; ...
      S0.eps_xy, S0.eps_yy, S0.eps_yz; ...
      S0.eps_xz, S0.eps_yz, S0.eps_zz ];

fprintf('Point (x,y,z) = (%.4f, %.4f, %.4f) m\n', x0,y0,z0);
fprintf('Strain tensor [eps_ij]:\n');
disp(E);

% ========= 2. 3x3 sub-image of the mid-thickness slice =========
N = 61;
xv = linspace(-params.Lx/2, params.Lx/2, N);
yv = linspace(-params.Ly/2, params.Ly/2, N);
zv = linspace(0, params.H, N);
[X,Y,Z] = meshgrid(xv,yv,zv);

S = strain_tensor(X,Y,Z,params);

% Find the layer closest to H/2
zmid = params.H/2;
[~,k] = min(abs(zv - zmid));

fields = { ...
    S.eps_xx(:,:,k), S.eps_xy(:,:,k), S.eps_xz(:,:,k); ...
    S.eps_xy(:,:,k), S.eps_yy(:,:,k), S.eps_yz(:,:,k); ...
    S.eps_xz(:,:,k), S.eps_yz(:,:,k), S.eps_zz(:,:,k) ...
    };

labels = { ...
    '\epsilon_{xx}','\epsilon_{xy}','\epsilon_{xz}'; ...
    '\epsilon_{yx}','\epsilon_{yy}','\epsilon_{yz}'; ...
    '\epsilon_{zx}','\epsilon_{zy}','\epsilon_{zz}' ...
    };

% Unified color scale range (for easy comparison)
all_abs = cellfun(@(A) max(abs(A(:))), fields);  % 3x3 array
cmax = max(all_abs(:));                          % <- Key: Take the maximum value of all elements

% Avoid illegal cmax (0, NaN, Inf)
if isempty(cmax) || ~isfinite(cmax) || cmax<=0
    cmax = 1e-12;
end


figure('Color','w','Position',[100 80 1200 1050]);
for i=1:3
    for j=1:3
        subplot(3,3,(i-1)*3+j);
        contourf(xv,yv,fields{i,j},40,'LineStyle','none');
        axis equal tight;
        xlabel('x (m)'); ylabel('y (m)');
        title(labels{i,j},'Interpreter','tex');
        colormap(parula); 
        colorbar;
        caxis([-cmax, cmax]);  % Unified color code
    end
end
sgtitle(sprintf('Strain tensor components on slice z=H/2=%.3f m', zmid));
