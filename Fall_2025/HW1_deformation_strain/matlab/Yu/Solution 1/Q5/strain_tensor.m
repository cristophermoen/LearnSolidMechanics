function S = strain_tensor(x, y, z, params)
% strain_tensor_prof: small-strain tensor from professor-based field
%

F = uvw_model(x,y,z,params);

% small strains (symmetric)
S.eps_xx = F.ux;
S.eps_yy = F.vy;
S.eps_zz = F.wz;

S.eps_xy = 0.5*(F.uy + F.vx);
S.eps_xz = 0.5*(F.uz + F.wx);
S.eps_yz = 0.5*(F.vz + F.wy);

% Engineering shear (gamma = 2*eps)
S.gamma_xy = 2*S.eps_xy;
S.gamma_xz = 2*S.eps_xz;
S.gamma_yz = 2*S.eps_yz;
end
