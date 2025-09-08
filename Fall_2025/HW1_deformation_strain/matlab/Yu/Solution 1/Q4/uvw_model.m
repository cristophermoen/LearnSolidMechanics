function out = uvw_model(x, y, z, params)
% uvw_prof_model: Extended professor-based displacement field 
%
% u(x,y,z) = A * sin(pi*y/H) * (2*x/Lx - 1)   
% v(x,y,z) = -(Delta/H) * y                   
% w(x,y,z) = A * sin(pi*y/H) * (2*z/Lz - 1)   
%
% Inputs:
%   x,y,z        : coordinates
%
% Output struct fields:
%   u,v,w        
%   ux,uy,uz ... 

Lx    = params.Lx;
H     = params.H;
Lz    = params.Lz;
Delta = params.Delta;
A     = params.A;

u = A .* sin(pi.*y./H) .* (2.*x./Lx - 1);   
v = -(Delta./H) .* y;                       
w = A .* sin(pi.*y./H) .* (2.*z./Lz - 1);   

% u gradients
ux = A .* sin(pi.*y./H) .* (2./Lx);
uy = A .* (pi./H) .* cos(pi.*y./H) .* (2.*x./Lx - 1);
uz = zeros(size(x));

% v gradients
vx = zeros(size(x));
vy = -(Delta./H) + 0.*x;
vz = zeros(size(x));

% w gradients
wx = zeros(size(x));
wy = A .* (pi./H) .* cos(pi.*y./H) .* (2.*z./Lz - 1);
wz = A .* sin(pi.*y./H) .* (2./Lz);

% --- output ---
out.u=u; out.v=v; out.w=w;
out.ux=ux; out.uy=uy; out.uz=uz;
out.vx=vx; out.vy=vy; out.vz=vz;
out.wx=wx; out.wy=wy; out.wz=wz;
end
