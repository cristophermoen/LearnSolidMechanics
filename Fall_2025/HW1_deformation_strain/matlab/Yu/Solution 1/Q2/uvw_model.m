function F = uvw_prof_model(x, y, z, params)
% uvw_prof_model (Professor-based, corner origin, y vertical)
% Displacement field + full 3x3 first derivatives
%
% u(x,y,z) = A*sin(pi*y/H)*(2x/Lx - 1)
% v(x,y,z) = -(Delta/H)*y
% w(x,y,z) = A*sin(pi*y/H)*(2z/Lz - 1)

Lx = params.Lx;   H = params.H;   Lz = params.Lz;
A  = params.A;    Delta = params.Delta;

% ----- Displacements -----
F.u = A .* sin(pi.*y./H) .* (2.*x./Lx - 1);
F.v = -(Delta./H) .* y;
F.w = A .* sin(pi.*y./H) .* (2.*z./Lz - 1);

% ----- First derivatives -----
% du/dx, du/dy, du/dz
F.ux = A .* sin(pi.*y./H) .* (2./Lx);
F.uy = A .* (pi./H) .* cos(pi.*y./H) .* (2.*x./Lx - 1);
F.uz = zeros(size(x));

% dv/dx, dv/dy, dv/dz
F.vx = zeros(size(x));
F.vy = -(Delta./H) .* ones(size(x));
F.vz = zeros(size(x));

% dw/dx, dw/dy, dw/dz
F.wx = zeros(size(x));
F.wy = A .* (pi./H) .* cos(pi.*y./H) .* (2.*z./Lz - 1);
F.wz = A .* sin(pi.*y./H) .* (2./Lz);
end
