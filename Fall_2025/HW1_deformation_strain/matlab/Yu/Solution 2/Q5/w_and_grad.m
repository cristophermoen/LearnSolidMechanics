function [w, dwdx, dwdy, dwdz] = w_and_grad(x, y, z, params)
% w_and_grad: Compute w(x,y,z) and its gradients for the bulging bearing model
%
% Inputs:
%   x,y,z   : coordinates (can be scalars or arrays of the same size)
%   params  : struct with fields：
%             .Lx, .Ly, .H   : bearing dimensions 
%             .Delta         : vertical compression at z=H 
%             .A             : bulging amplitude 
%
% Outputs:
%   w       : vertical displacement field 
%   dwdx    : ∂w/∂x
%   dwdy    : ∂w/∂y
%   dwdz    : ∂w/∂z

% ---- unpack parameters ----
Lx    = params.Lx;
Ly    = params.Ly;
H     = params.H;
Delta = params.Delta;
A     = params.A;

% ---- common factors ----
pi_over_H = pi / H;
pi_over_Lx = pi / Lx;
pi_over_Ly = pi / Ly;

% ---- basis functions ----
sin_z = sin(pi_over_H .* z);
cos_z = cos(pi_over_H .* z);
cos_x = cos(pi_over_Lx .* x);
sin_x = sin(pi_over_Lx .* x);
cos_y = cos(pi_over_Ly .* y);
sin_y = sin(pi_over_Ly .* y);

% ---- displacement ----
w = -(Delta/H) .* z + A .* sin_z .* cos_x .* cos_y;

% ---- gradients ----
dwdx = -A .* sin_z .* (pi_over_Lx) .* sin_x .* cos_y;
dwdy = -A .* sin_z .* (pi_over_Ly) .* cos_x .* sin_y;
dwdz = -(Delta/H) + A .* (pi_over_H) .* cos_z .* cos_x .* cos_y;
end
