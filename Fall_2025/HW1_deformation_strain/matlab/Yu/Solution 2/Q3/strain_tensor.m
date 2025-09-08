function S = strain_tensor(x,y,z,params)
% strain_tensor: Small-strain tensor from u,v,w fields
%
% Definitions:
%   eps_xx = du/dx,  eps_yy = dv/dy,  eps_zz = dw/dz
%   eps_xy = 0.5*(du/dy + dv/dx)
%   eps_xz = 0.5*(du/dz + dw/dx)
%   eps_yz = 0.5*(dv/dz + dw/dy)
% Engineering shear:
%   gamma_xy = 2*eps_xy,  etc.
%
% Output:
%   S.eps_xx, S.eps_yy, S.eps_zz, S.eps_xy, S.eps_xz, S.eps_yz
%   S.gamma_xy, S.gamma_xz, S.gamma_yz
%   S.Eeps : 3x3 tensor [eps_ij] 
%   S.Geng : 3x3 with engineering shear (off-diagonals = gamma/2 on Eeps)   
%
% Required：uv_model.m, w_and_grad.m

% --- lateral field u,v & grads ---
[u,v,dudx,dudy,dudz,dvdx,dvdy,dvdz] = uv_model(x,y,z,params);

% --- vertical field w & grads ---
[~, dwdx, dwdy, dwdz] = w_and_grad(x,y,z,params);

% --- small-strain components  ---
eps_xx = dudx;
eps_yy = dvdy;
eps_zz = dwdz;

eps_xy = 0.5*(dudy + dvdx);
eps_xz = 0.5*(dudz + dwdx);
eps_yz = 0.5*(dvdz + dwdy);

% --- engineering shear ---
gamma_xy = 2*eps_xy;
gamma_xz = 2*eps_xz;
gamma_yz = 2*eps_yz;

% Pack into struct 
S.eps_xx = eps_xx; S.eps_yy = eps_yy; S.eps_zz = eps_zz;
S.eps_xy = eps_xy; S.eps_xz = eps_xz; S.eps_yz = eps_yz;
S.gamma_xy = gamma_xy; S.gamma_xz = gamma_xz; S.gamma_yz = gamma_yz;

% Tensor matrix (symmetric)  3x3 
% [eps] = [ eps_xx  eps_xy  eps_xz
%           eps_xy  eps_yy  eps_yz
%           eps_xz  eps_yz  eps_zz ]
S.Eeps = cat(3, ...
    eps_xx, eps_xy, eps_xz, ...
    eps_xy, eps_yy, eps_yz, ...
    eps_xz, eps_yz, eps_zz);

% For convenience, also return a 3x3 "display" matrix per-point
S.Geng = S.Eeps; % same numeric values;
end
