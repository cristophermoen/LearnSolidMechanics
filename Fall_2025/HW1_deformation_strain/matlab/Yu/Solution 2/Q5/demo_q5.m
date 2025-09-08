clear; close all; clc;

% -------- Parameters --------
params.Lx = 0.10;    % m
params.Ly = 0.10;    % m
params.H  = 0.10;    % m

inch = 0.0254;
params.Delta = 2.0*inch;   % 2 inches compression
params.A  = 0.002;         % vertical bulging amplitude
params.B  = 0.001;         % lateral amplitude

% -------- Grid --------
N = 31; 
xv = linspace(-params.Lx/2, params.Lx/2, N);
yv = linspace(-params.Ly/2, params.Ly/2, N);
zv = linspace(0, params.H, N);
[X,Y,Z] = meshgrid(xv, yv, zv);

% -------- Strain tensor --------
S = strain_tensor(X,Y,Z,params);

% -------- Fields --------
fields = { ...
    S.eps_xx, S.eps_xy, S.eps_xz; ...
    S.eps_xy, S.eps_yy, S.eps_yz; ...
    S.eps_xz, S.eps_yz, S.eps_zz ...
    };

labels = { ...
    '\epsilon_{xx}','\epsilon_{xy}','\epsilon_{xz}'; ...
    '\epsilon_{yx}','\epsilon_{yy}','\epsilon_{yz}'; ...
    '\epsilon_{zx}','\epsilon_{zy}','\epsilon_{zz}' ...
    };

% ============================================================
%   fig.1：slice 
% ============================================================
figure('Color','w','Position',[50 50 1400 1200]);

for i=1:3
    for j=1:3
        subplot(3,3,(i-1)*3+j);
        val = fields{i,j};
        slice(X,Y,Z,val, ...
              [0], [0], [params.H/2]);
        shading interp;
        xlabel('x'); ylabel('y'); zlabel('z');
        title(labels{i,j},'Interpreter','tex');
        axis equal tight; view(3);
        colorbar;
    end
end
sgtitle('3D slice plots of strain tensor components');

% ============================================================
%   fig.2：isosurface
% ============================================================
figure('Color','w','Position',[100 100 1400 1200]);

for i=1:3
    for j=1:3
        subplot(3,3,(i-1)*3+j);
        val = fields{i,j};
        % Select the isosurface: take the median of the component
        midval = median(val(:));
        p = patch(isosurface(X,Y,Z,val,midval));
        isonormals(X,Y,Z,val,p);
        set(p,'FaceColor','red','EdgeColor','none','FaceAlpha',0.6);
        camlight headlight; lighting gouraud;
        xlabel('x'); ylabel('y'); zlabel('z');
        title([labels{i,j},' (iso)'],'Interpreter','tex');
        axis equal tight; view(3);
    end
end
sgtitle('3D isosurface plots of strain tensor components');
