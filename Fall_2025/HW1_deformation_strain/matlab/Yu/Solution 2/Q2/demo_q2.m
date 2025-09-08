% Parameters (consistent with Q1)
params.Lx = 0.20;   % m
params.Ly = 0.20;   % m
params.H  = 0.05;   % m
params.Delta = 0.010; % m
params.A  = 0.002;  % m

% Evaluate at center and mid-thickness
x0 = 0; y0 = 0; z0 = params.H/2;
[w0, dwdx0, dwdy0, dwdz0] = w_and_grad(x0, y0, z0, params)

% Construct a flat grid (z = H/2) to view the distribution (optional)
nx = 81; ny = 81;
xv = linspace(-params.Lx/2, params.Lx/2, nx);
yv = linspace(-params.Ly/2, params.Ly/2, ny);
[X,Y] = meshgrid(xv, yv);
Z = (params.H/2) * ones(size(X));
[W, Dx, Dy, Dz] = w_and_grad(X, Y, Z, params);

% Simple visualization of the distribution of w in the middle thickness (demo)
figure; surf(X, Y, W); xlabel('x'); ylabel('y'); zlabel('w at z=H/2');
title('w(x,y,z=H/2) with bulging');
shading interp; box on
