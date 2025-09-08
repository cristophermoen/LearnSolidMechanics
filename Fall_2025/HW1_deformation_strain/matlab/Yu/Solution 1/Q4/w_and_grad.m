function [w, wx, wy, wz] = w_and_grad(x, y, z, params)
% w_and_grad: Professor-based displacement field (corner origin)
%
% Inputs:
%   x,y,z       : coordinates (scalars/arrays of same size) 
%   params.H    : thickness along y
%   params.Lz   : strip length along z (m) 
%   params.A    : bulging amplitude (m) 
%
% Outputs:
%   w           : displacement field 
%   wx,wy,wz    : its gradients ∂w/∂x, ∂w/∂y, ∂w/∂z

H  = params.H;
Lz = params.Lz;
A  = params.A;

% Field
w  = A .* sin(pi.*y./H) .* (2.*z./Lz - 1);

% Gradients
wx = zeros(size(w));                                   % ∂w/∂x = 0
wy = A .* (pi./H) .* cos(pi.*y./H) .* (2.*z./Lz - 1);  % ∂w/∂y
wz = A .* sin(pi.*y./H) .* (2./Lz);                    % ∂w/∂z
end
