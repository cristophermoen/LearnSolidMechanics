function [u,v,dudx,dudy,dudz,dvdx,dvdy,dvdz] = uv_model(x,y,z,params)
% uv_model: Lateral (x,y) displacement field consistent with bulging mode
%
% Inputs:
%   params.B : lateral amplitude (m) 
Lx = params.Lx; Ly = params.Ly; H = params.H;
if isfield(params,'B'), B = params.B; else, B = 0; end

pi_over_H  = pi / H;
pi_over_Lx = pi / Lx;
pi_over_Ly = pi / Ly;

sin_z = sin(pi_over_H.*z);
cos_z = cos(pi_over_H.*z);
sin_x = sin(pi_over_Lx.*x);
cos_x = cos(pi_over_Lx.*x);
sin_y = sin(pi_over_Ly.*y);
cos_y = cos(pi_over_Ly.*y);

% Displacements 
u =  B .* sin_z .* sin_x .* cos_y;
v =  B .* sin_z .* cos_x .* sin_y;

% First derivatives 
dudx =  B .* sin_z .* (pi_over_Lx).*cos_x .* cos_y;
dudy = -B .* sin_z .* (pi_over_Ly).*sin_x .* sin_y;
dudz =  B .* cos_z .* (pi_over_H) .* sin_x .* cos_y;

dvdx = -B .* sin_z .* (pi_over_Lx).*sin_x .* sin_y;
dvdy =  B .* sin_z .* (pi_over_Ly).*cos_x .* cos_y;
dvdz =  B .* cos_z .* (pi_over_H) .* cos_x .* sin_y;
end
