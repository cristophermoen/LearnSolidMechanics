% demo_q3.m

clear; close all; clc;

%% ========= Parameters =========
params.Lx = 0.10;     
params.H  = 0.10;      
params.Lz = 0.10;    
inch = 0.0254;
params.Delta = 2.0*inch;   
params.A     = 0.002;      

%% ========= 1) Strain tensor at any point (user-defined point) =========

x0 = 0.05;                 % m
y0 = params.H/2;           % m
z0 = params.Lz/2;          % m

x0 = min(max(x0,0), params.Lx);
y0 = min(max(y0,0), params.H );
z0 = min(max(z0,0), params.Lz);

S0 = strain_tensor(x0, y0, z0, params);

E = [ S0.eps_xx, S0.eps_xy, S0.eps_xz; ...
      S0.eps_xy, S0.eps_yy, S0.eps_yz; ...
      S0.eps_xz, S0.eps_yz, S0.eps_zz ];

fprintf('Point (x,y,z) = (%.4f, %.4f, %.4f) m\n', x0, y0, z0);
fprintf('Strain tensor [eps_ij]:\n'); disp(E);
fprintf('Engineering shear (gamma):  gamma_xy = %.3e,  gamma_xz = %.3e,  gamma_yz = %.3e\n', ...
        S0.gamma_xy, S0.gamma_xz, S0.gamma_yz);

try
    lam = eig((E+E.')/2);   
    lam = sort(lam,'descend');
    fprintf('Principal strains (sorted):  [% .6e, % .6e, % .6e]\n', lam(1), lam(2), lam(3));
catch
end

%% ========= 2) 3×3 contour plot of the y = H/2 slice =========
Nx = 81; Nz = 81;                 
xv = linspace(0, params.Lx, Nx);
zv = linspace(0, params.Lz, Nz);
[X, Z] = meshgrid(xv, zv);
Y = params.H/2 * ones(size(X));   
S = strain_tensor(X, Y, Z, params);

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

figure('Color','w','Position',[80 80 1200 1050]);
for i=1:3
    for j=1:3
        Aij = fields{i,j};
        subplot(3,3,(i-1)*3+j);
        contourf(xv, zv, Aij', 40, 'LineStyle','none'); 
        axis equal tight; box on;
        xlabel('x (m)'); ylabel('z (m)');
        title(labels{i,j},'Interpreter','tex');
        cb = colorbar;
        cmax = max(abs(Aij(:)));
        if isempty(cmax) || ~isfinite(cmax) || cmax<=0, cmax=1e-12; end
        caxis([-cmax cmax]); cb.Ticks = linspace(-cmax,cmax,5);
    end
end
sgtitle(sprintf('Strain tensor components on plane y = H/2 = %.3f m', params.H/2));
